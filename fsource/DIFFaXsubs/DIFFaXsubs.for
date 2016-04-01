************************************************************************
*                                                                      *
*     Copyright 1987-2010 Michael M. J. Treacy and Michael W. Deem     *
*                                                                      *
************************************************************************
************************************************************************
*****************      Source file DIFFaXsubs.for       ****************
************************************************************************
************************************************************************
********************* version 1.813, 19th May, 2010 ********************
******* Trimmed for GSAS-II use by R. Von Dreele March, 2016 ***********
************************************************************************
************************************************************************
*
* ______________________________________________________________________
* Title: block data
* Author: MMJT
* Date: 23 Oct 1988
* Description: Sets up some global constants which help make the
* code easier to read.
*
*      COMMON VARIABLES:
*
*        modifies:  CENTRO, ELECTN, NEUTRN, NONE, GAUSS,
*                   LORENZ, PS_VGT, PV_GSS, PV_LRN, X_RAY
* ______________________________________________________________________
*
      block data
      implicit none
*
      integer*4 NONE, CENTRO, GAUSS, LORENZ, PS_VGT, PV_GSS, PV_LRN,
     |          X_RAY, NEUTRN, ELECTN
*
* The common block 'consts' also occurs in the file 'DIFFaX.inc'
      common /consts/ NONE, CENTRO, GAUSS, LORENZ, PS_VGT, PV_GSS,
     |                PV_LRN, X_RAY, NEUTRN, ELECTN
*
      data NONE, CENTRO /0, 1/
      data GAUSS, LORENZ, PS_VGT, PV_GSS, PV_LRN /1, 2, 3, 4, 5/
      data X_RAY, NEUTRN, ELECTN /0, 1, 2/
*
      end
*
* ______________________________________________________________________
* Title: AGLQ16
* Author: MWD
* Date: 18 Aug 1988
* Description:  This routine does adaptive 16-point Gauss-Legendre
* quadrature in an interval (h,k,a) to (h,k,b) of reciprocal space.
* ok is returned as .TRUE. if GLQ16 hasn't blown its stack. The
* integrated result is returned in AGLQ16.
*
*      ARGUMENTS:
*            h   -  reciprocal lattice vector h-component. (input).
*            k   -  reciprocal lattice vector k-component. (input).
*            a   -  l-value of the lower bound of reciprocal
*                   lattice integration region. (input).
*            b   -  l-value of the upper bound of reciprocal
*                   lattice integration region. (input).
*            ok  -  logical flag indicating all went well. (output).
*
*      AGLQ16 returns the adaptively integrated value.
* ______________________________________________________________________
*
      real*8 function AGLQ16(h, k, a, b, ok)
      include 'DIFFaX.par'
*     save
*
      integer*4 h, k
      real*8 a, b
      logical ok
*
      integer*4 maxstk, stp, n, n2
      parameter(maxstk = 200)
      real*8 sum, sum1, sum2, sum3, epsilon, epsilon2, GLQ16,
     |         stk(maxstk), d1, d2, d3, x
      parameter(epsilon = FIVE * eps4)
*
* external function
      external GLQ16
*
      AGLQ16 = ZERO
* initalize stack; top is at highest index
      stp = maxstk
* get first integration
      sum3 = GLQ16(h, k, a, b, ok)
      if(.not.ok) goto 999
      n2 = 2
      n = 0
      d1 = a
      d3 = b
      sum = ZERO
   20 d2 = HALF * (d1 + d3)
        n = n + 2
        sum1 = GLQ16(h, k, d1, d2, ok)
        if(.not.ok) goto 999
        sum2 = GLQ16(h, k, d2, d3, ok)
        if(.not.ok) goto 999
        x = sum1 + sum2
* determine figure of merit
        epsilon2 = max(epsilon, abs(x) * epsilon)
        if(abs(x - sum3).gt.epsilon2) then
* the area of these two panels is not accurately known
* check for stack overflow
          if(stp.lt.3) goto 980
* push right one onto stack
          stk(stp) = sum2
          stk(stp-1) = d2
          stk(stp-2) = d3
          stp = stp - 3
          d3 = d2
          sum3 = sum1
        else
* this panel has been accurately integrated; save its area
          sum = sum + x
* get next panel
* check for stack underflow--happens when no panels left to pop off
          if(stp.eq.maxstk) goto 30
          d3 = stk(stp+1)
          d1 = stk(stp+2)
          sum3 = stk(stp+3)
          stp = stp + 3
        endif
        if(n.eq.n2) then
          n2 = n2 * 2
        endif
      goto 20
   30 AGLQ16 = sum
      return
  980 write(op,402) 'Stack overflow in AGLQ16.'
      return
  999 write(op,402) 'GLQ16 returned an error to AGLQ16.'
      write(op,403) h, k, d1, d2
      return
  402 format(1x, a)
  403 format(3x,'at: h = ',i3,' k = ',i3,' l1 = ',g12.5,' l2 = ',g12.5)
      end
*
* ______________________________________________________________________
* Title: APPR_F
* Author: MMJT
* Date: 11 Feb 1989
* Description: This subroutine returns polynomial approximations of f
*  at 16 points, for use in GLQ16. The f's are returned in a MAX_L x 16
*  array. The order of the polynomial is n-1, where n is input. The 16
*  points in reciprocal space for which f's are needed are given by
*  (h,k,ag_l(i)), where h, k and ag_l are input.
*  The ll are the n sampling points, whose f's must be calculated
*  by direct calls to GET_F, and from which the interpolations are
*  calculated. list contains the indices of those n points. There is no
*  need to interpolate for these values since we know them already!
*
*      ARGUMENTS:
*            f      -  Array of layer form factors. (output).
*            h      -  reciprocal lattice vector h-component. (input).
*            k      -  reciprocal lattice vector k-component. (input).
*            ll     -  an array of n l-values to be used in the
*                      interpolation. (input).
*            ag_l   -  an array of 16 l_values at which interpolated
*                      form factors are required. (input).
*            n      -  the order of the polynomial approximation.
*                                                              (input).
*            list   -  A list of the indices of the n ll points
*                      entered. Interpolation is not needed at these
*                      input values. (input).
*            ok     -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:     a0, b0, c0, d0, n_layers
* ______________________________________________________________________
*
      subroutine APPR_F(f,h,k,ll,ag_l,n,list,ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 n, list(n)
      real*8 ll(n), ag_l(16)
      integer*4 h, k
      complex*16 f(MAX_L,16)
      logical ok
*
      logical know_f
      integer*4 i, j, m, p, max_poly
      parameter (max_poly = 10)
      real*8 Q2, l
      complex*16 ff(MAX_L,max_poly), fa(MAX_L), f_ans, f_error
*
* external subroutines (Some compilers need them declared external)
*      external POLINT, GET_F
*
* statement function
* Q2 is the value of 1/d**2 at hkl
      Q2(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
*
* sample GET_F n times
      do 10 i = 1, n
        call GET_F(ff(1,i), Q2(h,k,ll(i)), ll(i))
   10 continue
*
* store the n sampled f values for each layer i in fa.
* call POLINT for each layer in turn
* do this 16 times.
      do 20 m = 1, 16
* check to see that we haven't just calculated f at l(m)
        know_f = .false.
        do 30 i = 1, n
          if(m.eq.list(i)) then
            p = i
            know_f = .true.
          endif
   30   continue
* if we have, then simply copy it
        if(know_f) then
          do 40 i = 1, n_layers
            f(i,m) = ff(i,p)
   40     continue
        else
* else, use polynomial interpolation.
          do 50 i = 1, n_layers
            do 60 j = 1, n
              fa(j) = ff(i,j)
   60       continue
            call POLINT(ll,fa,n,ag_l(m),f_ans,f_error,ok)
            if(.not.ok) goto 999
            f(i,m) = f_ans
   50     continue
        endif
   20 continue
*
      return
  999 write(op,100) 'POLINT returned an error to APPR_F.'
      return
  100 format(1x, a)
      end
**
** ______________________________________________________________________
** Title: ATOMS
** Authors: MMJT
** Date: 18 Feb 90
** Description: This routine lists the legal atom names accepted by
** DIFFaX. These names are taken from the file 'sfname'.
**
**      ARGUMENTS:
**           No arguments are used.
**
**      COMMON VARIABLES:
**            uses:  sfname
**
**        modifies:  No COMMON variables are modified.
** ______________________________________________________________________
**
*      subroutine ATOMS()
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 i, LENGTH
*      character*4 atom_name
*      character*80 atomline
*      character*120 line
**
** external function
*      external LENGTH
**
** sfname has been opened once already,
** no need for elaborate error checking
*      open(unit = sf, file = sfname, status = 'old', err = 999)
**
*      write(op,200) 'Legal atom types are:'
**
** skip first few lines which contain no data, until we come to Hydrogen
*    1 read(sf,'(a)') line
*        if(line(1:4).ne.'H   ') goto 1
*      backspace(unit = sf, err = 500)
**
** initialize atomline
*      atomline = ' '
** write 10 atom names to a line
*   10 i = 0
*   20 read(sf, '(a)', end = 100) line
*        atom_name = line(1:4)
*        write(atomline(i*7+1:i*7+7),'(3a)') '''', atom_name, ''' '
*        i = i + 1
*        if(mod(i,10).eq.0) then
*          write(op,210) atomline
*          goto 10
*        endif
*        goto 20
*  100 continue
**
*      close(unit = sf, err = 600)
**
*  999 return
*  500 write(op,220) 'Unable to backspace scattering factor file ''',
*     |    sfname(1:LENGTH(sfname)), '''.'
*      return
*  600 write(op,220) 'Unable to close scattering factor file ''',
*     |    sfname(1:LENGTH(sfname)), '''.'
*      return
*  200 format(1x, a)
*  210 format(2x, a)
*  220 format(1x, 'ERROR in ATOMS: ', 3a)
*      end
**
* ______________________________________________________________________
* Title: BINPOW
* Author: MMJT
* Date: 18 Mar 1990
* Description:  This function breaks down the number 'n' into its
* binary representation. The result is stored in the global array 'pow'.
* This is used for efficiently multiplying a square matrix by itself
* n times if the RECURSIVE option was chosen for a finite number of
* layers.
*
* n must be such that n <= RCSV_MAX+1 <= 2**(MAX_BIN+1) - 1
*
*      ARGUMENTS:
*            n   -  number to binarize. (input).
*
*      COMMON VARIABLES:
*            uses:  max_pow
*
*        modifies:  max_pow, pow
*
* BINPOW returns logical TRUE if no problems were encountered.
* ______________________________________________________________________
*
      logical function BINPOW(n)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 n
*
      integer*4 i, j, itmp
*
      BINPOW = .false.
*
      itmp = n
      max_pow = 0
      i = MAX_BIN + 1
*
   10 i = i - 1
        j = itmp / 2**(i-1)
* j should be either 1 or zero if n was within bounds
        if(j.eq.1) then
          pow(i) = 1
          itmp = itmp - 2**(i-1)
          if(max_pow.eq.0) max_pow = i
        else if(j.eq.0) then
          pow(i) = 0
        else
          goto 999
        endif
      if(i.gt.1) goto 10
*
      BINPOW = .true.
      return
  999 write(op,100) 'ERROR in BINPOW: Invalid exponent ', n
      write(op,100) ' Maximum No. of layers allowed = ', RCSV_MAX
      write(op,100) ' Maximum supported by DIFFaX is', 2**(MAX_BIN+1)-1
      return
  100 format(1x, a, i8)
      end
*
* ______________________________________________________________________
* Title: BOUNDS
* Authors: MMJT
* Date: 24 Feb 1990
* Description: This function translates the value x so that it lies
* within the range 0 to 1.
*
*      ARGUMENTS:
*                  x  -  real number to convert. (input).
*
*      COMMON VARIABLES:
*
*                  No COMMON variables are used
*
*      BOUNDS returns the translated value of x. x is not modified.
* ______________________________________________________________________
*
      real*8 function BOUNDS(x)
      include 'DIFFaX.par'
*
      real*8 x
*
      real*8 y
*
      y = x - int(x) + ONE
      y = y - int(y)
*
* allow for rounding error
      if((ONE - y).lt.eps5) y = ZERO
*
      BOUNDS = y
*
      return
      end
*
* ______________________________________________________________________
* Title: CHK_SYM
* Author: MMJT
* Date: 15 Aug 1989; 21 Jan 1995
* Checks the user's assertions in the data file about the symmetry of
* the diffraction pattern. The symmetry is used to define the smallest
* possible volume in reciprocal space over which integration can
* occur and be representative of the whole pattern. CHK_SYM gives the
* user a crude 'goodness of fit', or tolerance, within which a random
* sampling of intensities fit the asserted symmetry. CHK_SYM does
* not override the user's judgment if the fit is poor, unless the cell
* dimensions or angle are inconsistent with the requested symmetry. The
* intention of CHK_SYM is to alert the user that the data does not
* conform to the expected symmetry. Useful for debugging datafiles.
*
*      ARGUMENTS:
*            ok  -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:  cell_a, cell_b, cell_gamma, max_angle, SymGrpNo,
*                   pnt_grp, tolerance, PI, PI2, PS_VGT, RAD2DEG
*
*        modifies:  max_var, h_mirror, k_mirror, check_sym, SymGrpNo
* ______________________________________________________________________
*
      subroutine CHK_SYM(ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical ok
*
      logical diad, triad, tetrad
      logical TST_ROT, TST_MIR
      logical cell90, cell120, eq_sides
      integer*4 GET_SYM, LENGTH, idum
      real*8 tmp
*
* external functions
      external TST_ROT, TST_MIR, GET_SYM, LENGTH
* external subroutine (Some compilers need them declared external)
*
* reinitialize random numbers in RAN3
      idum = -1
*
      diad   = .false.
      triad  = .false.
      tetrad = .false.
      cell90 =  abs(cell_gamma - HALF*PI) .lt. HALF*PI*eps6
      cell120 = abs(cell_gamma - PI2/THREE) .lt. PI2*eps6/THREE
*
* sample reciprocal space to get an idea of the sort of intensities
* that are out there.
* 360 degrees, no symmetry (-1)
* there is nothing to check. The center of symmetry is a given.
      if(SymGrpNo.eq.1) then
        max_var = ZERO
        goto 900
      endif
*
* 180 degrees, rotation only (2/M, 1st setting)
      if(SymGrpNo.eq.2) then
        diad = TST_ROT(2, idum, ok)
        if(.not.ok) goto 999
        goto 900
      endif
*
* 180 degrees, vertical mirror (2/M, 2nd setting)
      if(SymGrpNo.eq.3) then
        h_mirror = TST_MIR(1, idum, ok)
        if(.not.ok) goto 999
        tmp = max_var
        max_var = ZERO
        k_mirror = TST_MIR(2, idum, ok)
        if(.not.ok) goto 999
        max_var = max(tmp, max_var)
        goto 900
      endif
*
* 90 degrees, vertical mirrors (MMM)
      if(SymGrpNo.eq.4) then
        if(.not.cell90) goto 910
        diad = TST_ROT(2, idum, ok)
        if(.not.ok) goto 999
        tmp = max_var
        max_var = ZERO
        h_mirror = TST_MIR(1, idum, ok)
        if(.not.ok) goto 999
        tmp = max(tmp, max_var)
        max_var = ZERO
        k_mirror = TST_MIR(2, idum, ok)
        if(.not.ok) goto 999
        max_var = max(tmp, max_var)
        goto 900
      endif
*
* the following point groups require equi-sided cells
      eq_sides = abs(cell_a - cell_b) .le. HALF*eps6*(cell_a + cell_b)
      if(.not.eq_sides) goto 920
*
*
* 120 degrees, rotation only (-3)
      if(SymGrpNo.eq.5) then
        if(.not.cell120) goto 910
        triad = TST_ROT(3, idum, ok)
        if(.not.ok) goto 999
        tmp = max_var
        max_var = ZERO
        h_mirror = TST_MIR(1, idum, ok)
        if(.not.ok) goto 999
        tmp = max(tmp, max_var)
        max_var = ZERO
        hk_mirror = TST_MIR(3, idum, ok)
        if(.not.ok) goto 999
        max_var = max(tmp, max_var)
        goto 900
      endif
*
* 60 degrees, vertical mirrors (-3M)
      if(SymGrpNo.eq.6) then
        if(.not.cell120) goto 910
        triad = TST_ROT(3, idum, ok)
        if(.not.ok) goto 999
        tmp = max_var
        max_var = ZERO
        h_mirror = TST_MIR(1, idum, ok)
        if(.not.ok) goto 999
        tmp = max(tmp, max_var)
        max_var = ZERO
        hk_mirror = TST_MIR(3, idum, ok)
        if(.not.ok) goto 999
        max_var = max(tmp, max_var)
        goto 900
      endif
*
* 90 degrees, rotation (4/M)
      if(SymGrpNo.eq.7) then
        if(.not.cell90) goto 910
        tetrad = TST_ROT(4, idum, ok)
        if(.not.ok) goto 999
        goto 900
      endif
*
* 45 degrees, vertical mirrors (4/MMM)
      if(SymGrpNo.eq.8) then
        if(.not.cell90) goto 910
        tetrad = TST_ROT(4, idum, ok)
        if(.not.ok) goto 999
        tmp = max_var
        max_var = ZERO
        h_mirror = TST_MIR(1, idum, ok)
        if(.not.ok) goto 999
          if(.not.h_mirror) then
          tmp = max(tmp, max_var)
          max_var = ZERO
          hk_mirror = TST_MIR(3, idum, ok)
          if(.not.ok) goto 999
          endif
        max_var = max(tmp, max_var)
        goto 900
      endif
*
* 60 degrees, rotation (6/M)
      if(SymGrpNo.eq.9) then
        if(.not.cell120) goto 910
        diad = TST_ROT(2, idum, ok)
        if(.not.ok) goto 999
        tmp = max_var
        max_var = ZERO
        triad = TST_ROT(3, idum, ok)
        if(.not.ok) goto 999
        max_var = max(tmp, max_var)
        goto 900
      endif
*
* 30 degrees, vertical mirrors (6/MMM)
      if(SymGrpNo.eq.10) then
        if(.not.cell120) goto 910
        diad = TST_ROT(2, idum, ok)
        if(.not.ok) goto 999
        tmp = max_var
        max_var = ZERO
        triad = TST_ROT(3, idum, ok)
        if(.not.ok) goto 999
        tmp = max(tmp, max_var)
        max_var = ZERO
        h_mirror = TST_MIR(1, idum, ok)
        if(.not.ok) goto 999
        if(.not.h_mirror) then
          tmp = max(tmp, max_var)
          max_var = ZERO
          hk_mirror = TST_MIR(3, idum, ok)
          if(.not.ok) goto 999
        endif
        max_var = max(tmp, max_var)
      endif
*
      if(SymGrpNo.gt.10) goto 500
*
  900 write(op,202)
     | 'The diffraction data fits the point group symmetry ''',
     | pnt_grp(1:LENGTH(pnt_grp)),''''
      if(max_var.gt.eps6 .and. max_var.le.eps1) then
        write(op,203) '  with a tolerance of one part in ',
     |                  nint(ONE / max_var)
      else if(max_var.gt.eps1) then
        write(op,204) '  with a tolerance of one part in ',
     |                  ONE / max_var
      else
        write(op,100)
     |      '  with a tolerance better than one part in a million.'
      endif
*
  500 return
*
* The user's guess is inconsistent with cell_gamma.
* Override the user.
  910 write(op,200) 'The cell angle of',cell_gamma * RAD2DEG,
     |            ' degrees,'
      write(op,202) ' is inconsistent with point group symmetry ''',
     |            pnt_grp(1:LENGTH(pnt_grp)),''''
      write(op,300)
      check_sym = .false.
      SymGrpNo = GET_SYM(ok)
      if(.not.ok) goto 999
      write(op,205) pnt_grp(1:LENGTH(pnt_grp))
      if(tolerance.gt.eps6 .and. tolerance.le.eps1) then
        write(op,203) '  with a tolerance of one part in ',
     |                  nint(ONE / tolerance)
      else if(tolerance.gt.eps1) then
        write(op,204) '  with a tolerance of one part in ',
     |                  ONE / tolerance
      else
        write(op,100)
     |      '  with a tolerance better than one part in a million.'
      endif
      return
*
* The user's guess is inconsistent with cell dimensions.
* Override the user.
  920 write(op,201) 'The cell a and b dimensions, ',
     |            cell_a,' Angstroms by ',cell_b,' Angstroms,'
      write(op,202) '   are inconsistent with point group symmetry ''',
     |            pnt_grp(1:LENGTH(pnt_grp)),''''
      write(op,300)
* reset check_sym flag, since we are now evaluating from scratch
      check_sym = .false.
      max_var = ZERO
      SymGrpNo = GET_SYM(ok)
      if(.not.ok) goto 999
      write(op,205) pnt_grp(1:LENGTH(pnt_grp))
*
      if(tolerance.gt.eps6 .and. tolerance.le.eps1) then
        write(op,203) '  with a tolerance of one part in ',
     |                  nint(ONE / tolerance)
      else if(tolerance.gt.eps1) then
        write(op,204) '  with a tolerance of one part in ',
     |                  ONE / tolerance
      else
        write(op,100)
     |      '  with a tolerance better than one part in a million.'
      endif
*
      return
*
  999 write(op,100) 'ERROR in CHK_SYM'
      return
  100 format(1x, a)
  200 format(1x, a, f7.3, a)
  201 format(1x, 2(a, f7.3), a)
  202 format(1x, 3a)
  203 format(1x, a, i6)
  204 format(1x, a, f3.1)
  205 format(1x, 'Diffraction point symmetry is found to be ''',a,'''')
  300 format(1x, 'Re-evaluating diffraction symmetry')
      end
*
* ______________________________________________________________________
* Title: CHOICE
* Author: MWD
* Date: 18 Aug 1988 (modified by MMJT 6 August 1989)
* Description: This function compares a string with an array of keys,
* and returns the index of the key that matches the first characters
* of the string. CHOICE returns -1 if an error is found.
*
*      ARGUMENTS:
*            flag  -  Character string to be matched. (input).
*            list  -  Array of character strings. (input).
*            n     -  Length of the array 'list'. (input).
* ______________________________________________________________________
*
      integer*4 function CHOICE(flag, list, n)
      implicit none
*
      integer*4 n
*      character*(*) flag, list(n)
      character*(*) flag
      character*80 list(n)
*
      integer*4 i, j1, j2, LENGTH
*
* external function
      external LENGTH
*
      i = 1
      j1 = LENGTH(flag)
*
   10 j2 = LENGTH(list(i))
* see if the string contained in list(i) is identical to that in flag
      if(j1.eq.j2 .and. index(flag, list(i)(1:j2)).eq.1) goto 20
        i = i + 1
        if(i.le.n) goto 10
*
   20 if(i.gt.n) i = -1
      CHOICE = i
*
      return
      end
*
* ______________________________________________________________________
* Title: CHWDTH
* Author: MMJT
* Date: 6 Mar 1995; 31st Oct 1996
* Description:  This routine adds shape broadening caused by the finite
* lateral width of the layers. This routine does not add the broadening
* caused by finite width in the stacking direction. That is handled in
* a different manner by the routine INTEN2 and associated routines.
* The broadening handled here is added phenomenologically by applying a
* Lorentzian profile to the computed intensity at each peak in the h-k
* plane. If we are on the 00l axis the shape broadening is not
* symmetrical. For a Lorentzian shape broadening intensity profile in
* the h-k plane, the curvature of the Ewald sphere ensures a sharp
* onset of intensity with a slowly decaying tail at higher l values. For
* a symmetrical disk, the integrated intensity decays logarithmically.
* In the event that the crystal width is different in the h, and h-
* perpendicular directions, the disk of confusion becomes elongated
* into a streak perpendicular to l, and the tail becomes more
* Lorentzian-like. This is modelled in a phenomenological fashion, by
* mixing Lorentzian and logarithmic terms in a manner that depends
* on the ratio Wa/Wb. The costly logarithm function is avoided by
* using its derivative.
* When off the 00l axis, the broadening is modelled as a symmetric
* Lorentzian whose half-width depends on the angle that the Ewald 
* sphere intercepts the disk of confusion (controlled by l). If
* the lateral dimensions are not equal, then the half width is also
* dependent on h and k. The Lorentzian is pre-computed in OPTIMZ to gain
* computational speed.
* This routine is called by GETSPC.
*
*      ARGUMENTS:
*            h       -  Reciprocal space index. (input).
*            k       -  Reciprocal space index. (input).
*            l0      -  Lower bound of the l reciprocal space index
*                       that is being integrated over. (input).
*            l1      -  Lower bound of the l reciprocal space index
*                       that is being integrated over. (input).
*            x       -  The integrated intensity value along the
*                       line defined by h,k,l0,l1. (input).
*            m       -  The current index of the array 'spec'
*                       corresponding to h,k,l0,l1. (input).
*          max_indx  -  Maximum array value in spec that will be
*                       accessed. (input).
*
*      COMMON VARIABLES:
*            uses:      spec, brd_spec, FFACT_SIZE, formfactor, d_theta
*                       ffact_scale, ffhkcnst, ffwdth
*
*        modifies:      spec, brd_spec
* ______________________________________________________________________
*
      subroutine CHWDTH(h,k,l0,l1,x,m,max_indx)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k, m, max_indx
*
      real*8 l0, l1, x
*
      integer*4 n, p, i, indx
*
      real*8 S, h_wdth, n_hw, d_hk, norm, l, scale, avg, xx, dx, tmp
*
* indx indexes into the arrays spec and brd_spec
* n indexes into the array formfactor
* p is the index of the centroid of the array formfactor
* h_wdth contains the effective half-width of the size broadening
* ffhkcnst depends only on h and k, and was computed in GETSPC
* scale contains the calibration term for the formfactor array
* d_hk is the average radius in reciprocal Angstroms that the
* Ewald sphere of radius theta+d_theta intercepts the 00l plane
* brd_spc is used for temporary storage.
* We are only accessing half of the symmetrical formfactor array
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
*
************************************************************************
* special case if we are on the l axis
* Don't bother if the intensity is negligible in the first place
      if(h.eq.0 .and. k.eq.0 .and. x.gt.TEN*tiny_inty) then
        l = HALF*(l0+l1)
        d_hk = TWO*l*d_theta / (lambda*cell_c*ffwdth*ffwdth)
        norm = ZERO
        indx = 0
        xx = Wa / Wb
        if(xx.gt.ONE) xx = ONE / xx
   10   indx = indx + 1
          tmp = ONE / (ONE + indx*d_hk)
* balance streak, versus disk of confusion (phenomenological treatment)
* xx=1 means disk, xx=0 means streak. Intensity falls off more
* rapidly for a streak
          dx = ((ONE-xx)*sqrt(dble(indx))*tmp + xx)*tmp
          if(m+indx-1.le.max_indx) brd_spc(m+indx-1) = dx
          norm = norm + dx
* eps5 is reasonable. However, it may be worth experimenting more.    
          if(dx.lt.eps5) goto 20
          goto 10
   20   continue
*
      norm = x / norm
      do 30 i = 0, indx - 1
        if(m+i.le.max_indx) spec(m+i) = spec(m+i) + norm*brd_spc(m+i)
   30 continue
*
* We were on the 00l axis, we can exit now
        return
      endif
*
************************************************************************
* We are not on the l-axis. Broadening is handled differently.
* scale relates the formfactor array indices to the spec array indices.
* h_wdth and ffact_scale should never be zero. In case they are,
* make scale large enough so that n.ge.p-1 and the loop below is
* exited early
      h_wdth = ffhkcnst / sqrt(S(h,k,HALF*(l0+l1)))
      if(h_wdth.gt.ZERO .and. ffact_scale.gt.ZERO) then
        scale = d_theta / (ffact_scale * h_wdth)
      else
        scale = FFACT_SIZE
      endif
*
      p = FFACT_SIZE/2 + 1
      norm = ONE
      brd_spc(m) = ONE
      indx = 0
   40 indx = indx + 1
        n_hw = indx*scale
        n = int(n_hw)
        if(n.ge.p-1) goto 50
* linear interpolation of the pre-computed pseudo-Lorentzian
        xx = n_hw - n
        avg = (ONE-xx)*formfactor(p+n) + xx*formfactor(p+n+1)
        if(m+indx.le.max_indx) brd_spc(m+indx) = avg
        if(m-indx.gt.0)        brd_spc(m-indx) = avg
* intensity x is being redistributed. We will need to normalize later
        norm = norm + TWO*avg
        goto 40
   50 continue
*
      norm = x / norm
      spec(m) = spec(m) + norm*brd_spc(m)
      do 60 i = 1, indx - 1
        if(m+i.le.max_indx) spec(m+i) = spec(m+i) + norm*brd_spc(m+i)
        if(m-i.gt.0)        spec(m-i) = spec(m-i) + norm*brd_spc(m-i)
   60 continue
*
      return
      end
*
************************************************************************
****************************LINPACK ROUTINES****************************
************************************************************************
* 
* The following are the standard Linpack routines for solving complex
* simultaneous equations. They were found to reduce DIFFaX run time by
* significant amount (30% in one case) compared with the Numerical 
* Recipes routines LUDCMP and LUBKSB. The only changes are
*                      
*                         complex -> complex*16
*                         real    -> real*8
*                         real()  -> dble()
*                         aimag   -> dimag
*
************************************************************************
* ______________________________________________________________________
* Title: CGESL (LINPACK ROUTINE)
* Author: cleve moler, university of new mexico, argonne national lab.
* Date: linpack. this version dated 08/14/78
* Description:
*     CGESL solves the complex system
*     a * x = b  or  ctrans(a) * x = b
*     using the factors computed by cgeco or CGEFA.
*
*     on entry
*
*        a       complex(lda, n)
*                the output from cgeco or CGEFA.
*
*        lda     integer
*                the leading dimension of the array  a .
*
*        n       integer
*                the order of the matrix  a .
*
*        ipvt    integer(n)
*                the pivot vector from cgeco or CGEFA.
*
*        b       complex(n)
*                the right hand side vector.
*
*        job     integer
*                = 0         to solve  a*x = b ,
*                = nonzero   to solve  ctrans(a)*x = b  where
*                            ctrans(a)  is the conjugate transpose.
*
*     on return
*
*        b       the solution vector  x .
*
*     error condition
*
*        a division by zero will occur if the input factor contains a
*        zero on the diagonal.  technically this indicates singularity
*        but it is often caused by improper arguments or improper
*        setting of lda .  it will not occur if the subroutines are
*        called correctly and if cgeco has set rcond .gt. 0.0
*        or CGEFA has set info .eq. 0 .
*
*     to compute  inverse(a) * c  where  c  is a matrix
*     with  p  columns
*           call cgeco(a,lda,n,ipvt,rcond,z)
*           if (rcond is too small) go to ...
*           do 10 j = 1, p
*              call CGESL(a,lda,n,ipvt,c(1,j),0)
*        10 continue
*
*     subroutines and functions
*
*     blas CAXPY,CDOTC
*     fortran conjg
*
* ______________________________________________________________________
*
      subroutine CGESL(a,lda,n,ipvt,b,job)
      implicit none
*
      integer lda,n,ipvt(1),job
      complex*16 a(lda,1),b(1)
*
      complex*16 CDOTC,t
      integer k,kb,l,nm1
* MMJT: external subroutine
      external CAXPY
* external function
      external CDOTC
*
      nm1 = n - 1
      if (job .ne. 0) go to 50
*
*        job = 0 , solve  a * x = b
*        first solve  l*y = b
*
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call CAXPY(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
*
*        now solve  u*x = y
*
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call CAXPY(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
*
*        job = nonzero, solve  ctrans(a) * x = b
*        first solve  ctrans(u)*y = b
*
         do 60 k = 1, n
            t = CDOTC(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/conjg(a(k,k))
   60    continue
*
*        now solve ctrans(l)*x = y
*
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + CDOTC(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
* ______________________________________________________________________
* Title: CAXPY (LINPACK ROUTINE)
* Author: jack dongarra
* Date: linpack, 3/11/78
* Description: constant times a vector plus a vector.
* ______________________________________________________________________
*
      subroutine CAXPY(n,ca,cx,incx,cy,incy)
      implicit none
*
      complex*16 cx(1),cy(1),ca
      integer n,incx,incy
*
      integer i,ix,iy
*
      if(n.le.0)return
      if (abs(dble(ca)) + abs(dimag(ca)) .eq. 0.0 ) return
      if(incx.eq.1.and.incy.eq.1)go to 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
*
*        code for both increments equal to 1
*
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
* ______________________________________________________________________
* Title: CDOTC (LINPACK ROUTINE)
* Author: jack dongarra
* Date: linpack,  3/11/78.
* Description:
*     forms the dot product of two vectors, conjugating the first
*     vector.
* ______________________________________________________________________
*
      complex*16 function CDOTC(n,cx,incx,cy,incy)
      implicit none
*
      complex*16 cx(1),cy(1)
      integer incx,incy,n
*
      complex*16 ctemp
      integer i,ix,iy
*
      ctemp = (0.0,0.0)
      CDOTC = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = ctemp + conjg(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      CDOTC = ctemp
      return
*
*        code for both increments equal to 1
*
   20 do 30 i = 1,n
        ctemp = ctemp + conjg(cx(i))*cy(i)
   30 continue
      CDOTC = ctemp
      return
      end
* ______________________________________________________________________
* Title: CGEFA (LINPACK ROUTINE)
* Author: cleve moler, university of new mexico, argonne national lab.
* Date: linpack. this version dated 08/14/78
* Description:
*
*     CGEFA factors a complex matrix by gaussian elimination.
*
*     CGEFA is usually called by cgeco, but it can be called
*     directly with a saving in time if  rcond  is not needed.
*     (time for cgeco) = (1 + 9/n)*(time for CGEFA) .
*
*     on entry
*
*        a       complex(lda, n)
*                the matrix to be factored.
*
*        lda     integer
*                the leading dimension of the array  a .
*
*        n       integer
*                the order of the matrix  a .
*
*     on return
*
*        a       an upper triangular matrix and the multipliers
*                which were used to obtain it.
*                the factorization can be written  a = l*u  where
*                l  is a product of permutation and unit lower
*                triangular matrices and  u  is upper triangular.
*
*        ipvt    integer(n)
*                an integer vector of pivot indices.
*
*        info    integer
*                = 0  normal value.
*                = k  if  u(k,k) .eq. 0.0 .  this is not an error
*                     condition for this subroutine, but it does
*                     indicate that CGESL or cgedi will divide by zero
*                     if called.  use  rcond  in cgeco for a reliable
*                     indication of singularity.
*
*     subroutines and functions
*
*     blas CAXPY,CSCAL,ICAMAX
*     fortran abs,aimag,real
* ______________________________________________________________________
*
      subroutine CGEFA(a,lda,n,ipvt,info)
      implicit none
*
      integer lda,n,ipvt(1),info
      complex*16 a(lda,1)
*
      complex*16 t
      integer ICAMAX,j,k,kp1,l,nm1
*
      complex*16 zdum
      real*8 cabs1
*
* MMJT: external subroutine
      external CSCAL, CAXPY
*
* external function
      external ICAMAX
* statement function
      cabs1(zdum) = abs(dble(zdum)) + abs(dimag(zdum))
*
*     gaussian elimination with partial pivoting
*
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
*
*        find l = pivot index
*
         l = ICAMAX(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
*
*        zero pivot implies this column already triangularized
*
         if (cabs1(a(l,k)) .eq. 0.0e0) go to 40
*
*           interchange if necessary
*
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
*
*           compute multipliers
*
            t = -(1.0e0,0.0e0)/a(k,k)
            call CSCAL(n-k,t,a(k+1,k),1)
*
*           row elimination with column indexing
*
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call CAXPY(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0e0) info = n
      return
      end
* ______________________________________________________________________
* Title: CSCAL (LINPACK ROUTINE)
* Author: jack dongarra
* Date: linpack,  3/11/78.
* Description: scales a vector by a constant.
* ______________________________________________________________________
*
      subroutine  CSCAL(n,ca,cx,incx)
      implicit none
*
      complex*16 ca,cx(1)
      integer incx,n
*
      integer i,nincx
*
      if(n.le.0)return
      if(incx.eq.1)go to 20
*
*        code for increment not equal to 1
*
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
*
*        code for increment equal to 1
*
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
* ______________________________________________________________________
* Title: ICAMAX (LINPACK ROUTINE)
* Author: jack dongarra
* Date: linpack, 3/11/78
* Description:
*     finds the index of element having max. absolute value.
* ______________________________________________________________________
*
      integer function ICAMAX(n,cx,incx)
      implicit none
*
      complex*16 cx(1)
      integer incx,n
*
      real*8 smax
      integer i,ix
      complex*16 zdum
*
* statement function
      real*8 cabs1
      cabs1(zdum) = abs(dble(zdum)) + abs(dimag(zdum))
*
      ICAMAX = 0
      if( n .lt. 1 ) return
      ICAMAX = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
*
*        code for increment not equal to 1
*
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         ICAMAX = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
*
*        code for increment equal to 1
*
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         ICAMAX = i
         smax = cabs1(cx(i))
   30 continue
      return
      end
*
* ______________________________________________________________________
* Title: CNTARG
* Author: MMJT
* Date: 24 Feb 1990
* Description: Counts the number of arguments in the character
* string 'line', and returns the answer in CNTARG. Legal separators
* are, blanks, tabs, and commas.
*
*      ARGUMENTS:
*            line  -  Line of characters to be parsed. (input).
*
*      CNTARG returns the number of integer arguments in the line.
*      Returns -1 if an error occurred.
* ______________________________________________________________________
*
      integer*4 function CNTARG(line)
      implicit none
*
      character*(*) line
*
      logical in_arg
      integer*4 blank, tab, comma, i, j, arg_cnt, lin_len
      parameter (tab=9, blank=32, comma=44)
*
      in_arg = .false.
      arg_cnt = 0
      CNTARG = -1
*
      lin_len = len(line)
      do 10 i = 1, lin_len
        j = ichar(line(i:i))
        if(j.eq.tab .or. j.eq.blank .or. j.eq.comma) then
          in_arg = .false.
        else
          if(.not.in_arg) then
            in_arg = .true.
            arg_cnt = arg_cnt + 1
          endif
        endif
   10 continue
*
* There cannot be more than (lin_len+1)/2 arguments
      if(arg_cnt.le.((lin_len+1)/2)) CNTARG = arg_cnt
*
      return
      end
*
* ______________________________________________________________________
* Title: DETUN
* Author: MMJT
* Date: 27 Jan 1989
* Description: This subroutine detunes the sharp peaks so that
* they can be integrated. In effect, it modifies the stacking
* probability matrix such that the determinant will never be
* exactly zero.
*
*      ARGUMENTS:
*            No input arguments.
*
*      COMMON VARIABLES:
*            uses:    n_layers
*
*        modifies:    detune
* ______________________________________________________________________
*
      subroutine DETUN()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 i, j
      real*8 delta
*
      delta = eps3
*
* A value of delta = 0.001 is found to be optimum.
* If preferred, user can specify 'delta' interactively
* using the following three lines.
C   30 write(op,400) 'Enter detune parameter'
C      read(cntrl,*,err=30,end=999) delta
C      if(CFile) write(op,401) delta
*
      do 10 i = 1, n_layers
        do 20 j = 1, n_layers
          detune(j,i) = ONE - abs(delta)
   20   continue
   10 continue
*
      return
C  999 stop 'ERROR: Bad delta value. DIFFaX aborted.'
C  400 format(1x, a)
C  401 format(1x, g12.5)
      end
*
* ______________________________________________________________________
** Title: DUMP
** Author: MWD and MMJT
** Date: 18 Aug 1988; 15 Mar 1995
** Description: This subroutine prints out the data read in from
** the data file. It reassures the user that the data was read in
** correctly.
**
**      ARGUMENTS:
**            infile  -  The name of the input data file. (input).
**            ok      -  logical flag indicating all went well.
**                                                      (output).
**
**      COMMON VARIABLES:
**            uses:      CFile, GAUSS, LORENZ, NEUTRN, PI2, PS_VGT,
**                       RAD2DEG, SymGrpNo, X_RAY, a_B, a_name,
**                       a_number, a_occup, a_pos, a_type, blurring,
**                       cell_a, cell_b, cell_c, cell_gamma, cntrl
**                       d_theta, e_sf, FWHM, inf_thick, l_actual
**                       l_alpha, l_cnt, l_g, l_n_atoms, l_r, l_symmetry
**                       lambda, l_seq, n_actual, n_atoms, n_layers
**                       pnt_grp, pv_gamma, pv_u, pv_v, pv_w, r_B11
**                       r_B12, r_B22, r_B23, r_B31, r_B33, rad_type
**                       recrsv, rndm, th2_max, th2_min
**                       tolerance, xplcit
**
**        modifies:      no COMMON variables are modified
** ______________________________________________________________________
**
*      subroutine DUMP(infile,ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      character*(*) infile
*      logical ok
**
*      real*8 scale, atom_cnt(MAX_TA), cum_atom_cnt(MAX_TA), norm
*      integer*4 i, i2, j, n, type, num_types, tot_types
*      integer*4 LENGTH, print_width
*      character*31 dmpfile
*      character*80 list(5)
**
** external function
*      external LENGTH
** external subroutine (Some compilers need them declared external)
**      external GETFNM
**
*      call GETFNM(infile, dmpfile, 'dmp', ok)
*      if(.not.ok) then
*        write(op,200) 'DUMP aborted'
*        goto 999
*      endif
*      if(dp.ne.op) open(unit = dp , file = dmpfile, status = 'new')
**
*      i2 = 0
*   99 write(op,200) 'Enter 1 for full atomic position dump'
*      read(cntrl,*,err=99,end=9999) i2
*      if(CFile) write(op,'(1x,i3)') i2
*      write(op,300) 'Writing data dump to file ''',
*     |               dmpfile(1:LENGTH(dmpfile)),'''. . .'
** sundry details about layers
*      write(dp,125) 'Number of layers = ', n_layers
*      write(dp,125) 'Number of unique layers = ', n_actual
*      write(dp,125) 'Number of different atom types = ', n_atoms
** cell dimensions
*      write(dp,140)
*     |    ' cell_a,      cell_b,      cell_c,      cell_gamma ',
*     |      cell_a, cell_b, cell_c, RAD2DEG * cell_gamma
** in-plane layer widths
*      if(finite_width) then
*        if(Wa.lt.inf_width) then
*          write(dp,126) 'width along a', Wa
*          write(dp,128) Wa/cell_a
*        else
*          write(dp,127) 'a'
*        endif
*        if(Wb.lt.inf_width) then
*          write(dp,126) 'width along b', Wb
*          write(dp,128) Wb/cell_b
*        else
*          write(dp,127) 'b'
*        endif
*      else
*        write(dp,200)
*     |    'Layers are to be treated as having infinite lateral width.'
*      endif
** radiation type
*      list(1) = 'X-RAY '
*      list(2) = 'NEUTRON '
*      list(3) = 'ELECTRON '
*      write(dp,100) 'radiation type = ', list(rad_type+1)
** wavelength
*      write(dp,170) 'Wavelength lambda = ', lambda
** symmetry
*      write(dp,100)
*     |'Diffraction point group symmetry specified = ',pnt_grp
*      if(SymGrpNo.eq.UNKNOWN) write(dp,175)
*     |'Tolerance on intensities in symmetry evaluation = +/-',
*     |                    tolerance * HUNDRED, ' %'
** instrumental broadening to simulate
*      list(1) = 'NONE'
*      list(2) = 'GAUSSIAN'
*      list(3) = 'LORENTZIAN'
*      list(4) = 'PSEUDO-VOIGT'
*      i = blurring + 1
** see if it's a pure Gaussian pseudo-Voigt
*      if(blurring.eq.PV_GSS) i = 2
** see if it's a pure Lorentzian pseudo-Voigt
*      if(blurring.eq.PV_LRN) i = 3
*      write(dp,100)
*     |  'Instrumental broadening to be simulated = ', list(i)
*      if(blurring.eq.GAUSS .or. blurring.eq.LORENZ) write(dp,170)
*     |'Full width half-maximum = ', FWHM
*      if(blurring.eq.PS_VGT) then
*        write(dp,200)
*     |'Pseudo-Voigt parameters:   u,       v,       w,     gamma'
*        write(dp,180) pv_u, pv_v, pv_w, pv_gamma
*      endif
*      if(blurring.eq.PV_GSS) then
*        write(dp,200) 'Gaussian parameters:      u,       v,       w'
*        write(dp,185) pv_u, pv_v, pv_w
*      endif
*      if(blurring.eq.PV_LRN) then
*        write(dp,200) 'Lorentzian parameters:    u,       v,       w'
*        write(dp,185) pv_u, pv_v, pv_w
*      endif
** are we to trim out the origin?
*      if(blurring.ne.NONE) then
*        if(trim_origin) then
*          write(dp,200) '  Intensity near origin will be ignored'
*        else
*          write(dp,200) '  Intensity near origin will be included'
*        endif
*      endif
**
** equivalent layers
*      write(dp,200) ' '
*      do 10 i = 1, n_layers
*        write(dp,110) 'LAYER ',i,
*     |            ' is equivalent to fundamental LAYER ',l_actual(i)
*        write(dp,170) '   Existence probability = ', l_g(i)
*   10 continue
**
*        do 11 j = 1, MAX_TA
*          cum_atom_cnt(j) = ZERO
*   11   continue
**
** layer structure
*      write(dp,200) ' '
*      do 30 i = 1, n_actual
*        write(dp,130) 'fundamental LAYER ',i
*        list(1) = 'NONE '
*        list(2) = 'CENTROSYMMETRIC '
*        write(dp,100) 'symmetry = ', list(l_symmetry(i)+1)
** write out the layer composition
*        write(dp,200) 'LAYER composition is:'
*        do 14 j = 1, MAX_TA
*          atom_cnt(j) = ZERO
*   14   continue
*        num_types = ZERO
*        scale = ONE
*        if(l_symmetry(i).eq.CENTRO) scale = TWO
*        do 15 j = 1, l_n_atoms(i)
*          type = a_type(j,i)
*          if(type.gt.num_types) num_types = type
*          atom_cnt(type) = atom_cnt(type) + a_occup(j,i) * scale
*   15   continue
*        if(num_types.gt.tot_types) tot_types = num_types
** accumulate weighted atom count
*        do 16 j = 1, MAX_TA
*          cum_atom_cnt(j) = cum_atom_cnt(j) + atom_cnt(j) * l_g(i)
*   16   continue
*        do 17 j = 1, num_types
*          write(dp,280) atom_l(j), ':', atom_cnt(j)
*   17   continue
*        write(dp,200) ' '
**
** do we want all the details about each atom?
*        if(i2.ne.1) goto 30
** yes, go for it.
*        do 20 j = 1, l_n_atoms(i)
*            write(dp,200) ' '
*          write(dp,120) 'ATOM number', a_number(j,i), ' in layer',i,
*     |     ' is ''', a_name(j,i),''''
*          write(dp,121)
*          write(dp,122) a_pos(1,j,i), a_pos(2,j,i), a_pos(3,j,i),
*     |      a_B(j,i), a_occup(j,i)
*          if(rad_type.eq.X_RAY) then
*            write(dp,140) 'X-ray scattering factor data Ai, Bi',
*     |                         (x_sf(n,a_type(j,i)),n=1,9)
*          else if(rad_type.eq.NEUTRN) then
*            write(dp,140) 'Neutron scattering factor data',
*     |                         n_sf(a_type(j,i))
*          else if(rad_type.eq.ELECTN) then
*            write(dp,140)'Electron scattering factor data Ai, Bi and Z',
*     |               (x_sf(n,a_type(j,i)),n=1,9), e_sf(a_type(j,i))
*          else
*            write(op,200) 'ERROR: Illegal radiation type in DUMP'
*          endif
*   20   continue
*        write(dp,200) ' '
*   30 continue
**
** Print out average composition
*      if(.not.xplcit) then
*        norm = ZERO
*        do 25 j = 1, num_types
*          norm = norm + cum_atom_cnt(j)
*   25   continue
*        write(dp,200) ' '
*        write(dp,200) 'Average crystal composition is:'
*        do 26 j = 1, num_types
*          write(dp,281) atom_l(j), ':', cum_atom_cnt(j) / norm
*   26   continue
*        write(dp,200) ' '
*      endif
**
** stacking details
*      write(dp,200) ' '
*      write(dp,200) 'STACKING'
*      if(recrsv) then
*        write(dp,210) 'RECURSIVELY'
*        if(inf_thick) then
*          write(dp,220) 'INFINITE'
*        else
*          write(dp,230) l_cnt
*        endif
*      else if(xplcit) then
*        if(rndm) then
*          write(dp,240)
*        else
*          write(dp,250)
*        endif
*        write(dp,260) 'Sequence for ', l_cnt, ' layers is:'
*        print_width = 25
*        j = 1
*        n = print_width
*   35   if(n.gt.l_cnt) n = l_cnt
*        write(dp,270) (l_seq(i), i = j, n)
*        if(n.lt.l_cnt) then
*          j = j + print_width
*          n = n + print_width
*          goto 35
*        endif
*      endif
**
** stacking transition data
*      write(dp,200) ' '
*      write(dp,200) ' '
*      write(dp,200) 'TRANSITIONS'
*      write(dp,200)
*     |      'Layer stacking probabilities and stacking vectors.'
*      write(dp,200)
*     |'  i  j |   alpha_ij      Rx_ij        Ry_ij        Rz_ij'
*      write(dp,200)
*     |'  -----|-------------------------------------------------'
*      do 40 i = 1, n_layers
*        write(dp,200)
*     |'       |'
*        do 50 j = 1, n_layers
*          write(dp,150) i, j, l_alpha(j,i),
*     |                   l_r(1,j,i), l_r(2,j,i), l_r(3,j,i)
*   50   continue
*   40 continue
**
*      write(dp,200) ' '
*      write(dp,200)
*     |      'Anisotropic Layer stacking ''uncertainty'' factors.'
** MMJT: 3/15/95. Ordering of B23 and B31 swapped. B31 -> C13
*      write(dp,200)
*     |'  i  j |     C11      C22      C33      C12      C13      C23'
*      write(dp,200)
*     |'  -----|-------------------------------------------------------'
*      do 60 i = 1, n_layers
*        write(dp,200)
*     |'       |'
*        do 70 j = 1, n_layers
** MMJT: 3/15/95. Ordering of B23 and B31 swapped
*          write(dp,151) i, j, r_B11(j,i), r_B22(j,i), r_B33(j,i),
*     |                        r_B12(j,i), r_B31(j,i), r_B23(j,i)
*   70   continue
*   60 continue
**
*      if(dp.ne.op) close(unit = dp)
*  999 return
* 9999 ok = .false.
*      write(op,200) 'DUMP aborted'
*      return
*  100 format(1x, 2a)
*  110 format(1x, a, i4, a, i4)
*  120 format(1x, a, i3, a, i2, 3a)
*  121 format(4x,'x_rel',8x,'y_rel',8x,'z_rel',10x,'DW',10x,'Occ')
*  122 format(1x, 5g13.5)
*  125 format(1x, a40, i4)
*  126 format(1x, 'layer characteristic ', a, ' = ', f9.2, ' Angstroms')
*  127 format(1x, 'layer characteristic ', a, ' = INFINITY Angstroms')
*  128 format(1x, '   which is equivalent to ', f9.2, ' unit cells')
*  130 format(1x, a, i4)
*  140 format(1x, a, / 5g13.5, / 4g13.5, i4)
*  150 format(1x, 2i3, ' |', 4(1x, g12.5))
*  151 format(1x, 2i3, ' |', 6(1x, f8.3))
*  170 format(1x, a, g12.5)
*  175 format(1x, a, g12.5, a)
*  180 format(21x, 4(2x, f7.3))
*  185 format(21x, 3(2x, f7.3))
*  200 format(1x, a)
*  210 format(1x, 'Stacking is to be treated ', a, ' by DIFFaX.')
*  220 format(1x, 'Number of layers along the fault direction is ', a)
*  230 format(1x, 'There are ',i5 ,' layers along the fault direction')
*  240 format(1x, 'Sequencing is defined RANDOMLY by DIFFaX.')
*  250 format(1x, 'Sequencing is defined EXPLICITLY by the user.')
*  260 format(1x, a, i4, a)
*  270 format(1x, 30i3)
*  280 format(23x, 2a, f7.2)
*  281 format(23x, 2a, f7.4)
*  300 format(1x, 3a)
*      end
*
* ______________________________________________________________________
* Title: EQUALB
* Author: MMJT
* Date: 13 Mar 1990; 21 July 1997; 3 July 2005
* Description:  This routine determines if all of the stacking
* uncertainty parameters are identical. There are six arrays to be
* tested, namely r_B11, r_B22, r_B33, r_B12, r_B23 and r_B31. These are
* passed one at a time from OPTIMZ as r_B. The average value of the
* r_B parameters is returned in a_B.
* The test fails if the user changes the sign, but not the amplitude, of
* some of the B11, B22 or B33. EQUALB returns false, and the calculation
* is then done the hard way.
*
*      ARGUMENTS:
*            r_B  -  Array of stacking uncertainty parameters. (input).
*            av_B  -  Average of r_B. (output).
*
*      COMMON VARIABLES:
*            uses the array 'there'. n_layers
* ______________________________________________________________________
*
      logical function EQUALB(r_B, av_B)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 r_B(MAX_L,MAX_L)
      real*8 av_B
*
      integer*4 i, j, m
      real*8 error
*
      av_B = ZERO
      m = 0
      do 10 i = 1, n_layers
        do 20 j = 1, n_layers
* Examine only those transitions that actually occur
          if(there(j,i)) then
            m = m + 1
            av_B = av_B + r_B(j,i)
          endif
   20   continue
   10 continue
*
* Take average
      if(m.ne.0) av_B = av_B / dble(m)
*
      error = ZERO
* find absolute deviation of coefficients
      do 30 i = 1, n_layers
        do 40 j = 1, n_layers
          if(there(j,i)) error = error + abs(r_B(j,i) - av_B)
   40   continue
   30 continue
*
      EQUALB = abs(error).le.abs(eps3*av_B)
*
      return
      end
*
* ______________________________________________________________________
** Title: FNDCTL
** Author: MMJT
** Date: 6 June 1989
** Description: Checks whether or not there is a control file, and
** sets the i/o unit specifier 'cntrl' to 99.
**
**      ARGUMENTS:
**           No arguments are used. All data is in 'COMMON'.
**
**      COMMON VARIABLES:
**            uses:   cfname, CFile
**
**        modifies:   cntrl
**
**      FNDCTL returns logical .true. if a control file was found in the
**      current directory.
** ______________________________________________________________________
**
*      logical function FNDCTL()
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 io_err, LENGTH
**
** external function
*      external LENGTH
**
** See if control file is present
*      inquire(file = cfname, exist = CFile)
** get data from ControlFile
*      if(CFile) then
*        cntrl = 99
*        open(unit = cntrl, file = cfname, status = 'old',
*     |      err = 999, iostat = io_err)
*      else
** get data from user entries to the screen
*        cntrl = ip
*      endif
*      FNDCTL = .true.
**
*      return
*  999 write(op,100)
*     |'A control file named ''', cfname(1:LENGTH(cfname)), ''''
*      write(op,101) 'exists, but there were problems opening it.'
*      write(op,102) 'iostat = ', io_err
*      FNDCTL = .false.
*      return
*  100 format(1x, 3a)
*  101 format(1x, a)
*  102 format(1x, a, i5)
*      end
**
** ______________________________________________________________________
** Title: GETFIL
** Author: MMJT
** Date: 5 June 1989
** Description: Checks for the presence of the structure data file and
** the structure factor data file in the current directory.
**
**      ARGUMENTS:
**            fname   -  The name of the input data file. (output).
**            ok      -  logical flag indicating all went well.
**                                                      (output).
**            ending  -  returns logical .true. if 'end' was read
**                       from the control file. (output).
**
**      COMMON VARIABLES:
**            uses:       cntrl, sfname, CFile
**
**        modifies:       no COMMON variables are modified
** ______________________________________________________________________
**
*      subroutine GETFIL(fname, ok, ending)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      character*(*) fname
*      logical ok, ending
**
*      character*31 tmp_nam
*      character*200 line
*      character*80 list(1)
*      integer*4 unit_no, LENGTH, CHOICE
**
*      integer*4 i, io_err
**
** external function
*      external LENGTH, CHOICE
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR
**
*      unit_no = cntrl
** get file name from ControlFile
*      write(op,104) 'Enter name of structure file.'
** use GETLNE to strip off any leading blanks etc...
*      call GETLNE(unit_no, line, *993)
*      write(fname,'(a)') line(1:len(fname))
*      if(CFile) write(op,104) fname
** check that the string does not say 'end'. First convert to upper case
** we must be sure that the string 'end' is not part of a legitimate
** file name. Use CHOICE for this.
*      tmp_nam = fname
*      call TOUPPR(fname)
*      list(1) = 'END '
*      i = CHOICE(fname, list, 1)
*      if(i.eq.1) then
*        ending = .true.
*        goto 999
*      endif
*      fname = tmp_nam
**
** Open structure data file.
*      write(op,100) 'Looking for structure data file ''',
*     |            fname(1:LENGTH(fname)),''''
*      open(unit = df, file = fname, status = 'old',
*     |            err = 991, iostat = io_err)
*      write(op,100) 'Opening structure file ''',
*     |       fname(1:LENGTH(fname)),''''
**
**Open scattering factor data file.
*      write(op,100)
*     | 'Looking for scattering factor data file ''',
*     |            sfname(1:LENGTH(sfname)),''''
*      open(unit = sf, file = sfname, status = 'old',
*     |            err = 992, iostat = io_err)
*      write(op,100) 'Opening scattering factor data file ''',
*     |       sfname(1:LENGTH(sfname)),''''
**
*      goto 999
*  991 write(op,100)
*     |'Problems opening structure data file ''',
*     |            fname(1:LENGTH(fname)),''''
*      write(op,103) 'IOSTAT = ',io_err
*      ok = .false.
*      goto 999
*  992 write(op,100)
*     |'Problems opening scattering factor file ''',
*     |            sfname(1:LENGTH(sfname)),''''
*      write(op,103) 'IOSTAT = ',io_err
*      ok = .false.
*      goto 999
*  993 write(op,104) 'ERROR reading filename in GETFIL'
*      ok = .false.
*  999 return
*  100 format(1x, 3a)
*  103 format(1x, a, i5)
*  104 format(1x, a)
*      end
**
** ______________________________________________________________________
** Title: GETFNM
** Author: MMJT
** Date: 13 Oct 1988
** Description: This subroutine makes sure that we do not write over
** an existing file, and creates a suitable filename (name2). name2
** equals 'name1.append' append is supplied by the calling routine.
** If 'name1.append' is the name of an existing file, then append is
** set to 'append#' where #=1,2.
** CLIP is the maximum length of the appendage. Thus if append
** = 'fred', and CLIP is 3, then the appendages become 'fre',
** 'fr1', 'fr2', ...'f10', 'f11' ... up to 'f99'.
** If CLIP is greater than 10, there is no restriction.
** NOTE: It is up to the user to ensure that filename lengths are
** compatible with the operating system that DIFFaX is being run on.
** 'CLIP' and 'MAX_NAM' are set in 'DIFFaX.par'
**
**      ARGUMENTS:
**            name1   -  The name of the input data file. (input).
**            name2   -  A derivative filename. (output).
**            append  -  A token to be appended to name2. (input).
**            ok      -  logical flag signalling all went well.
**                                                      (output).
** ______________________________________________________________________
**
*      subroutine GETFNM(name1,name2,append,ok)
*      include 'DIFFaX.par'
**     save
**
*      character*(*) name1, name2, append
*      logical ok
**
*      logical fexist
*      character*31 appnd2
*      integer*4 applen, i, myclip
*      integer*4 LENGTH
**
** external function
*      external LENGTH
** external subroutine (Some compilers need them declared external)
**      external NAMER
**
*      ok = .true.
*      fexist = .true.
*      appnd2 = append
*      applen = LENGTH(append)
*      myclip = CLIP
*      if(myclip.le.0) then
*        write(op,302) 'Illegal CLIP value ', CLIP
*        write(op,300) 'Default CLIP value of 3 is being assumed'
*        myclip = 3
*      endif
**
*      i = 0
*   10 call NAMER(name1,name2,appnd2)
*      inquire(file = name2, exist = fexist)
*      if(fexist) then
*        i = i + 1
*        if(i.lt.10) then
*          if(myclip.gt.10) then
*            write(appnd2(applen+1:applen+1),200) i
*          else
*            write(appnd2(myclip:myclip),200) i
*          endif
*        else if(i.ge.10 .and. i.lt.100) then
*          if(myclip.gt.10) then
*            write(appnd2(applen+1:applen+2),201) i
*          else if(myclip.gt.1) then
*            write(appnd2(myclip-1:myclip),201) i
*          else
*            ok = .false.
*            goto 999
*          endif
*        else if(i.ge.100 .and. i.lt.1000) then
*          if(myclip.gt.10) then
*            write(appnd2(applen+1:applen+3),202) i
*          else if(myclip.gt.2) then
*            write(appnd2(myclip-2:myclip),202) i
*          else
*            ok = .false.
*            goto 999
*          endif
*        else
*          ok = .false.
*          goto 999
*        endif
*        goto 10
*      endif
**
*      return
*  999 write(op,300) 'No new files can be created. . .'
*      write(op,300) '      . . . some old files should be deleted.'
*      write(op,301)
*     |      'Last file created was ''',name2(1:LENGTH(name2)),''''
*      return
*  200 format(i1)
*  201 format(i2)
*  202 format(i3)
*  300 format(1x, a)
*  301 format(1x, 3a)
*  302 format(1x, a, i5)
*      end
**
** ______________________________________________________________________
* Title: GETLAY
* Author: MMJT
* Date: 4 Oct 1989
* Description: This function generates a random sequence of layers
* weighted by the stacking probabilities. Needed only when the
* 'EXPLICIT' and 'RANDOM' options were specified in the 'STACKING'
* description.
*
*      ARGUMENTS:
*           No arguments are used. All data is in 'COMMON'.
*
*      COMMON VARIABLES:
*            uses:      l_cnt, l_g, n_layers, l_seq, l_alpha
*
*        modifies:      l_seq
*
*      GETLAY returns logical .true. if all went well.
* ______________________________________________________________________
*
      logical function GETLAY()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical okay
      character*80 messge
      integer*4 i, j, idum
      real*8 RAN3, x, sum
* external function
      external RAN3
*
      GETLAY = .false.
      okay = .true.
*
* initialize random numbers in RAN3
      idum = -1
*
      write(op,100)
     |      'Generating a random sequence of ', l_cnt, ' layers.'
*
* Get the first layer. Even though some stacking transition
* probabilities to and from the layer are non-zero, the layer itself
* may not actually exist anywhere in the crystal! Must use the
* existence probabilities, l_g.
      x = RAN3(idum)
      if(x.eq.ONE) x = x - eps7
      sum = ZERO
      i = 1
   10 sum = sum + l_g(i)
      if(x.gt.sum) then
        i = i + 1
        if(i.le.n_layers) goto 10
* if all goes well, we should not reach here.
        messge = 'GETLAY could not generate the first layer.$'
        goto 999
      endif
      l_seq(1) = i
*
* Now generate the remaining l_cnt-1 layers
      j = 2
   20 x = RAN3(idum)
      if(x.eq.ONE) x = x - eps7
      sum = ZERO
      i = 1
   30 sum = sum + l_alpha(i,l_seq(j-1))
      if(x.gt.sum) then
        i = i + 1
        if(i.le.n_layers) goto 30
* if all goes well, we should not reach here.
        write(messge,101) 'GETLAY could not generate layer ', j, '.$'
        goto 999
      endif
      l_seq(j) = i
      j = j + 1
      if(j.le.l_cnt) goto 20
*
      GETLAY = .true.
      return
  999 write(op,102) messge(1:index(messge,'$')-1)
      return
  100 format(1x, a, i4, a)
  101 format(a, i4, a)
  102 format(1x, 'ERROR: ', a)
      end
*
* ______________________________________________________________________
** Title: GETLNE
** Author: MWD  & MMJT
** Date: MWD: Original version 18 Aug 1988; MMJT: This version 9 Apr 1992
** Description: This function returns a non-blank, left-justified line,
** stripped of any comments, from unit 'unit_no'. Any error generates
** a return to the passed linenumber, presumably an error-handling
** routine in the calling procedure. Multiple comments, and nested
** braces on one line can be handled.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that data is to
**                        be read from. (input).
**            line     -  Line of data read as a character string.
**                                                          (output).
**            *        -  The label of a line in the calling routine
**                        to return to if an error occurs. (input).
** ______________________________________________________________________
**
*      subroutine GETLNE(unit_no, line, *)
*      include 'DIFFaX.par'
**     save
**
*      integer*4 unit_no
*      character*(*) line
**
*      character*1 tab, blank
** nominal line length is 200. Leave generous room for overflow.
*      character*250 rawline
*      integer*4 i, j, n, offset, lin_len
*      integer*4 l_br, r_br, first_left, last_right
**
*      tab = char(9)
*      blank = char(32)
*      lin_len = len(line)
**
*    1 read(unit_no, '(a)', err = 200, end = 200) rawline
*      line = rawline
**
** First check if user has gone beyond lin_len characters
*      do 2 i = len(rawline), lin_len + 1, -1
*        if(rawline(i:i).ne.blank) goto 998
*    2 continue
**
** strip out all tabs. Replace with blanks.
*    3 n = index(line,tab)
*      if(n.eq.0) goto 4
*        write(line(n:n),'(a)') blank
*        goto 3
*    4 continue
**
** Search for matching braces. Keep looping until all braces are found
*   10   first_left = 0
*        last_right = 0
*        l_br = 0
*        r_br = 0
*        i = 0
*   20   i = i + 1
*          if(i.gt.lin_len) goto 30
*          if(line(i:i).eq.'{') then
*            if(first_left.eq.0) first_left = i
*            l_br = l_br + 1
*          else if(line(i:i).eq.'}') then
*            last_right = i
*            r_br = r_br + 1
*          endif
*          if(r_br.gt.l_br) goto 999
*          if(l_br.eq.0 .or. l_br.gt.r_br) goto 20
*   30   continue
*        if(r_br.ne.l_br) goto 999
**
** if no more comments, then left justify the line
*        if(l_br.eq.0) goto 50
** else, remove comments, and shunt data to the left
*        offset = last_right - first_left
*        i = first_left
*        do 60 j = last_right + 1, lin_len
*          line(i:i) = line(j:j)
*          i = i + 1
*   60   continue
** blank out the end of the modified string
*        do 70 i = lin_len-offset, lin_len
*          line(i:i) = blank
*   70   continue
** repeat until no more braces are found
*      goto 10
**
** left adjust line
*   50 do 80 i = 1, lin_len
*        if((line(i:i).ne.blank)) then
*          do 90 j = 1, lin_len-i+1
*            line(j:j) = line(i+j-1:i+j-1)
*   90     continue
*          do 100 j = lin_len-i+2,lin_len
*            line(j:j) = blank
*  100     continue
*          goto 110
*        endif
*   80 continue
** If, after all this, the line is empty, repeat with next next line. 
*  110 if(line(1:1).eq.blank) goto 1
**
*      return
**
*  998 write(op,401) 'ERROR: Data line exceeds ',lin_len, ' characters.'
*      return 1
*  999 write(op,400) 'ERROR: Unbalanced braces in data file.'
*  200 return 1
**
*  400 format(1x, a)
*  401 format(1x, a, i3, a)
*      end
**
** ______________________________________________________________________
* Title: GETSAD
* Author: MMJT
* Date: 29 Oct 1989; 16th April 1999
* Description: This subroutine generates the selected area diffraction
* pattern (SADP). GETSAD returns .true. if all went well. The pattern
* is stored in the linear array 'spec'. 'spec' is read by WRTSAD
* to write the diffraction pattern image to disk as a binary file.
*
*      ARGUMENTS:
*            FN      -  Function name passed by reference. The
*                       choice is between GLQ16 (non-adaptive
*                       Gauss-Legendre integration), and AGLQ16
*                       (adaptive Gauss-Legendre integration). (input).
*            view    -  Choice of beam direction (input).
*                              1  =  normal to the plane k = 0
*                              2  =  normal to the plane h = 0
*                              3  =  normal to the plane h = k
*                              4  =  normal to the plane h = -k
*            l_upper -  Upper limit of l. (input).
*            hk_lim  -  Upper limit of h (or k). (output).
*            infile  -  The name of the input data file. (input).
*            ok      -  logical flag indicating all went well.
*                                                      (output).
*
*      COMMON VARIABLES:
*            uses:      a0, b0, c0, d0, lambda, has_l_mirror, sadblock,
*                       spec, loglin, brightness, X_RAY, rad_type
*
*        modifies:      scaleint
* ______________________________________________________________________
*
      subroutine GETSAD(FN, view, l_upper, hk_lim, infile, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 FN, l_upper
      integer*4 view, hk_lim
      character*(*) infile
      logical ok
*
      integer*4 h, k, i, j, n, info_step, info, cnt, LENGTH, origin
      real*8 x, S, S_value, ANGLE, W4, PNTINT, theta, Q2, l
      real*8 l_lower, dl, high1, high2, intervals
      parameter (intervals = TWENTY)
*
* external functions (FN is either GLQ16 or AGLQ16)
      external FN, LENGTH, PNTINT
*
* external subroutines (Some compilers need them declared external)
*      external XYPHSE, PRE_MAT
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
* ANGLE is the Bragg angle (in radians) of the h,k,l plane
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
* W4 is the X-ray polarization factor
      W4(theta) = HALF * (ONE + (cos(TWO*theta))**2)
*
      Q2 = FOUR / (lambda**2)
*
* check angles are meaningful
      S_value = S(0,0,l_upper)
      if(S_value.le.ZERO) then
        write(op,101)
     |     'ERROR: Illegal value in GETSAD: 1/d**2 = ',S_value
        ok = .false.
        goto 990
      endif
      if(S_value.gt.Q2) then
        write(op,100)
        l_upper = TWO / (lambda*sqrt(c0)) - eps10
        write(op,101) 'Upper bound reduced to ', l_upper
        S_value = S(0,0,l_upper)
      endif
*
* Set h and k limit
      if(view.eq.1) then
        hk_lim =  int(l_upper * sqrt(c0 / a0))
      else if(view.eq.2) then
        hk_lim =  int(l_upper * sqrt(c0 / b0))
      else if(view.eq.3) then
        hk_lim =  int(l_upper * sqrt(c0 / (a0 + b0 + d0)))
      else if(view.eq.4) then
        hk_lim =  int(l_upper * sqrt(c0 / (a0 + b0 - d0)))
      endif
*
* Get l increment
      dl = l_upper / dble(SADSIZE/2)
* reset l_upper so that our integration window, of width dl,
* straddles l = 0.
      l_upper = l_upper - HALF*dl
*
* Get diffraction pattern.
* First check scan limits, and make sure we are not going to overflow
* the array 'spec' which will be used to hold the scan information
* prior to writing to file in WRTSAD.
      if(has_l_mirror) then
        l_lower = ZERO - HALF*dl
        sadblock = SADSIZE/2
        if((hk_lim + 1)*(SADSIZE/2) .gt. MAX_SP) then
          hk_lim = MAX_SP/(SADSIZE/2) - 1
        endif
      else
        l_lower = -l_upper
        sadblock = SADSIZE
        if((hk_lim + 1)*SADSIZE .gt. MAX_SP) then
          hk_lim = MAX_SP/SADSIZE - 1
        endif
      endif
      info_step = nint( l_upper / (dl * intervals) )
      if(info_step.le.0) info_step = 1
*
      cnt = 0
      do 10 i = 0, hk_lim
* set h and k
        if(view.eq.1) then
          h =  i
          k =  0
        else if(view.eq.2) then
          h =  0
          k =  i
        else if(view.eq.3) then
          h =  i
          k =  i
        else if(view.eq.4) then
          h =  i
          k = -i
        else
          write(op,105) 'ERROR: Illegal view value in GETSAD', view
          goto 990
        endif
*
* write progress to screen
        write(op,102) h, k, infile(1:LENGTH(infile))
*
        call XYPHSE(h, k)
        call PRE_MAT(h, k)
*
* The following piece of monkey business is here because if there is
* no mirror then the l-loop only calculates SADSIZE-1 values. If there
* is a mirror then the l-loop makes SADSIZE/2 calculations. We wish
* the final binary data file to be in a square array SADSIZE x SADSIZE,
* so we pad out the first value with a zero.
*
        if(.not.has_l_mirror) then
          cnt = cnt + 1
          spec(cnt) = ZERO
        endif
*
        info = 0
        do 10 l = l_lower, l_upper-eps10, dl
* If we are at the origin, don't integrate, we'll probably overflow.
          if(i.eq.0 .and. abs(l+dl).le.dl+eps10) then
            x = ZERO
            origin = cnt + 1
          else if(S(h,k,l+dl).gt.Q2) then
            x = ZERO
          else
            x = FN(h,k,l,l+dl,ok)
            if(.not.ok) goto 999
            if(rad_type.eq.X_RAY) x = x*W4(ANGLE(h,k,l+HALF*dl))
          endif
          cnt = cnt + 1
* make sure we do not overflow
          if(cnt.gt.MAX_SP) goto 998
          spec(cnt) = x
          if(mod(info,info_step).eq.0) then
            if(loglin.eq.0) then
              if(ONE+x.gt.ZERO) then
* write out log(1+x), since x may get to be quite small
                write(op,103) 'l = ',l,' log(intensity) = ',log(ONE+x)
              else
                write(op,103) 'l = ', l, ' log(intensity) = ', ZERO
              endif
            else
              write(op,103) 'l = ', l, ' intensity = ', x
            endif
          endif
          info = info + 1
   10 continue
* check cnt
      if(cnt.lt.2) goto 980
*
* patch the intensity at the origin to put a well-defined peak there
      if(has_l_mirror) then
* origin = 1, make origin slightly bigger than proceeding value
        spec(origin) = (ONE+eps4) * spec(origin+1)
      else
* make it slightly larger than the biggest adjacent value
        spec(origin) = (ONE+eps4) * max(spec(origin-1),spec(origin+1))
      endif
*
* we need to find the second highest peak intensity so as to scale the
* data. The highest intensity should be at the origin, which we discard.
      high1 = ZERO
      high2 = ZERO
      do 30 i = 0, hk_lim
        do 30 j = 1, sadblock-1
          n = i*sadblock + j
* check scaling type. If logarithmic, make sure values are always +ve
          if(loglin.eq.0) then
            if(ONE+spec(n).gt.ZERO) then
              spec(n) = log(ONE + spec(n))
            else
              spec(n) = ZERO
            endif
          endif
* check if origin is the first value. If so, preset high value.
          if(n.eq.1 .and. origin.eq.1) then
            high1 = spec(origin)
            goto 30
          endif
          x = spec(n)
          if(j.eq.1) then
            if(x.gt.spec(n+1)) then
              if(x.gt.high1) then
                high2 = high1
                high1 = x
              else if(x.gt.high2) then
                high2 = x
              endif
            endif
          else
            if(x.gt.spec(n-1) .and. x.gt.spec(n+1)) then
              if(x.gt.high1) then
                high2 = high1
                high1 = x
              else if(x.gt.high2) then
                high2 = x
              endif
            endif
          endif
   30 continue
*
      if(loglin.ne.0 .and. high2.le.ZERO) then
        write(op,101)
     |  'ERROR in intensity scaling in GETSAD. ''scale factor'' = ',
     |                                                      high2
        ok = .false.
        goto 990
      endif
*
      if(loglin.eq.0 .and. high1.le.ZERO) then
        write(op,101)
     |  'ERROR in intensity scaling in GETSAD. ''scale factor'' = ',
     |                                                      high1
        ok = .false.
        goto 990
      endif
*
* If logarithmic, scale to the brightest peak
* If linear, scale to the 2nd brightest peak
* Note: Intensity scaling can be modified by the user-defined
* brightness value
      if(loglin.eq.0) then
        scaleint = brightness * (maxsad - ONE) / high1
      else
        scaleint = brightness * (maxsad - ONE) / high2
      endif
*
  990 return
  980 write(op,105)
     |     'Error in GETSAD: loop counter is too small. cnt = ', cnt
      ok = .false.
      return
  998 write(op,104) 'ERROR in GETSAD: spectrum array overflow at h = ',
     |                                    h,', k = ',k,', l = ',l
      ok = .false.
      return
  999 write(op,104) 'ERROR in GETSAD at h = ',h,', k = ',k,', l = ',l
      return
  100 format(1x, 'Upper bound exceeds 180 degrees!')
  101 format(1x, a, g12.5)
  102 format(1x, 'h = ', i3, ' k = ', i3, 10x, '''', a, '''')
  103 format(1x, a, f10.5, a, g12.5)
  104 format(1x, 2(a, i3), a, f10.5)
  105 format(1x, a, i3)
      end
*
* ______________________________________________________________________
* Title: GETSPC
* Authors: MWD and MMJT
* Date: 17 Mar 1989; Last tweaked on 7 Mar 1995
* Description: This subroutine calculates the spectrum.
*
*      ARGUMENTS:
*            FN      -  Function name passed by reference. The
*                       choice is between GLQ16 (non-adaptive
*                       Gauss-Legendre), and AGLQ16
*                       (adaptive Gauss-Legendre). (input).
*            infile  -  The name of the input data file. (input).
*
*      COMMON VARIABLES:
*            uses:      CFile, ELECTN, NEUTRN, SymGrpNo, X_RAY, a0
*                       any_sharp, b0, bnds_wt, c0, cntrl, d0, d_theta
*                       full_brd, lambda, mltplcty, rad_type, rot_only,
*                       spec, th2_max, th2_min, theta1, theta2
*
*        modifies:      full_shrp
*
*      GETSPC returns logical .true. if all went well.
* ______________________________________________________________________
*
      logical function GETSPC(FN,infile)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 FN
      character*(*) infile
*
      logical ok, SHARP, on_bndry, l_axis, shrp
      integer*4 h, k, h_lower, h_upper, k_lower, k_upper
      integer*4 m, i, max_indx
      integer*4 LENGTH
      real*8 S, Q, theta, tmp, tmp2, tmp3, fact, h_val, k_val
      real*8 HKANGL, LL, ANGLE, AGLQ16
      real*8 l, hk_th, x, GLQ16, l_max, min_th, max_th
      real*8 W1, l1, l0, d_l, INTENS, L_STEP, W2, W3, l00
      complex*16 f(MAX_L)
*
* external functions
      external FN, GLQ16, AGLQ16, INTENS, SHARP, L_STEP, LENGTH
* external subroutines (Some compilers need them declared external)
*      external XYPHSE, PRE_MAT, GET_F, CHWDTH
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
* LL is the maximum allowable l value for a given h,k and theta
      LL(theta,h,k) = sqrt((fact * (sin(theta))**2
     |                    - h*h*a0 - k*k*b0 - h*k*d0)/ c0)
* ANGLE is the Bragg angle (in radians) of the h,k,l plane
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
* HKANGL is the angle between the vector (h_val,k_val,0) and (1,0,0)
      HKANGL(k_val,h_val) = atan2(k_val*sqrt(a0*b0 - d0*d0*QUARTER),
     |                              (h_val*a0 + k_val*d0*HALF))
* These factors are for powder diffraction patterns
* W1 is the polarization factor for weighting x-ray intensities
* it also incorporates the Lorentz factor
      W1(theta) = (ONE + cos(TWO*theta) * cos(TWO*theta)) /
     |                  (sin(theta) * sin(TWO*theta))
* W2 is the neutron weighting factor--the Lorentz factor
      W2(theta) = ONE / (sin(theta) * sin(TWO*theta))
* W3 is the electron weighting factor--the Lorentz factor
      W3(theta) = ONE / (sin(theta) * sin(TWO*theta))
*
      GETSPC = .FALSE.
      ok = .true.
*
* Make sure we are within bounds. If not, adjust.
      min_th = HALF * th2_min
      max_th = HALF * th2_max
      max_indx = int((max_th - min_th) / d_theta + 1)
      if(max_indx.gt.MAX_SP) then
        d_theta = (max_th - min_th) / (MAX_SP - 1)
        max_indx = int((max_th - min_th) / d_theta + 1)
        write(op,300) ''
        write(op,250) 'd_theta is too small and has been adjusted to ',
     |                   TWO*d_theta*RAD2DEG
      endif
*
 1234 write(op,300)
     | 'Enter 1 for adaptive quadrature over all l values'
      write(op,300) 'on rows with "sharp" spots'
      read(cntrl,*,err=1234) full_shrp
      if(CFile) write(op,400) full_shrp
* zero out spectra
      do 10 i = 1, MAX_SP
        spec(i) = ZERO
   10 continue
* See if there is a chance of any sharp peaks.
* If so, an appropriate step along l is found, and any_sharp is .true.
      d_l = L_STEP(ok)
      if(.not.ok) goto 999
* If no sharp peaks were unambiguously detected, override user.
      if(d_l.eq.ZERO) full_shrp = 1
* determine extreme values of h
      Q = TWO * sin(max_th) / lambda
      fact = TWO / lambda
      fact = fact * fact
* h_upper is the largest value that the index h can have,
* consistent with the value of Q. (When cell_gamma is not 90 degrees,
* k is not necessarily zero at this extreme).
* In case the incredible happens, immune system to rescue.
      tmp3 = FOUR * a0 * b0 - d0 * d0
      if(tmp3.le.ZERO) goto 990
      tmp3 = TWO * Q * sqrt(ONE / tmp3)
      h_upper = int(tmp3 * sqrt(b0))
      h_lower = -h_upper
      k_upper = int(tmp3 * sqrt(a0))
      k_lower = -k_upper
* scan along h-axis from h_lower to h_upper
      do 20 h = h_lower, h_upper
* determine limits along k for a given h
        do 30 k = k_lower, k_upper
* if out of bounds, cycle
          if(S(h,k,ZERO).gt.Q*Q) goto 30
          l_axis = h.eq.0 .and. k.eq.0
          hk_th = theta1
          if(.not.l_axis) hk_th = HKANGL(dble(k), dble(h))
* see if in wedge to be scanned
          if((theta1-hk_th)*(theta2-hk_th).le.eps3 .or.
     |            SymGrpNo.eq.1) then
* if rotational symmetry only, do not take points on upper wedge plane
            if(rot_only .and. (theta2-hk_th).le.eps3 .and.
     |            SymGrpNo.ne.1) goto 30
            if(SymGrpNo.eq.11 .and. .not.l_axis) goto 30
            write(op,200) 'Integrating along l at ',h,k,
     |            '''',infile(1:LENGTH(infile)),''''
            on_bndry = abs(hk_th-theta1).le.eps3 .or.
     |                 abs(hk_th-theta2).le.eps3
* set up the phases in the structure factors and stacking vectors
* which depend on h and k only
            call XYPHSE(h, k)
            call PRE_MAT(h, k)
* assign a corrected shape-broadening half-width
            if(finite_width) then
              tmp2 = (h + k*cos(PI-cell_gamma))/(Wa*sin(PI-cell_gamma))
              tmp3 = k / Wb
              ffhkcnst = QUARTER*lambda*sqrt(a0*tmp2*tmp2+b0*tmp3*tmp3)
            endif
* get starting value of l
            if(l_axis) then
              tmp = min(d_theta, max_th)
              if(tmp.lt.min_th) tmp = min_th
              l1 = LL(tmp, h, k)
              shrp = any_sharp
            else
              tmp = ANGLE(h, k, ZERO)
              if(tmp.lt.min_th) then
                l1 = LL(min_th,h,k)
                tmp = ANGLE(h, k, l1)
              else
                l1 = ZERO
              endif
              if(any_sharp .and. full_shrp.ne.1) then
                shrp = SHARP(h, k, d_l)
                if(.not.ok) goto 999
              else
                shrp = any_sharp
              endif
            endif
* m indexes into the array spec
            m = int((tmp - min_th) / d_theta) + 1
            if(.not.shrp .or. full_shrp.eq.1) then
* broad streak or full adaptive integration over sharp spots
              if(full_shrp.eq.1 .or. full_brd.eq.1)
     |                write(op,300) 'Full adaptive integration'
* integrate each d_theta's range of reciprocal space
              do 40 theta = tmp, max_th-eps10, d_theta
                l0 = l1
                tmp2 = min(d_theta, max_th-theta)
                l1 = LL(theta+tmp2, h, k)
* sharp spots; do not use knowledge of where they are
                if(shrp) then
                  x = AGLQ16(h, k, l0, l1, ok)
                else
* broad streaks
                  x = FN(h, k, l0, l1, ok)
                endif
                if(.not.ok) goto 110
*
* include weighting factors for radiation type
                if(rad_type.eq.X_RAY) then
                  x = TWO * x * W1(theta + HALF * tmp2)
                else if(rad_type.eq.NEUTRN) then
                  x = TWO * x * W2(theta + HALF * tmp2)
                else if(rad_type.eq.ELECTN) then
                  x = TWO * x * W3(theta + HALF * tmp2)
                else
                  ok = .false.
                  goto 130
                endif
*
* see if not on l-axis
                if(.not.l_axis) then
* apply multiplicity factor
                  x = x * mltplcty
* if on boundary, apply appropriate weighting (mirror vs rotation only)
                  if(on_bndry) x = x * bnds_wt
                endif
                if(finite_width) then
                  call CHWDTH(h,k,l0,l1,x,m,max_indx)
                else
                  spec(m) = spec(m) + x
                endif
                m = m + 1
   40         continue
            else
* line of sharp spots--detuned delta functions
* use knowledge of where spots are
* make sure we do all l values a multiple of d_l
* starting with spot on hk-plane
              write(op,300) 'which is a line of sharp spots'
              l00 = ZERO
              if(l_axis) then
                l00 = d_l
   50           if(l00.ge.th2_min) goto 60
                  l00 = l00 + d_l
                  goto 50
   60           continue
              endif
              l_max = LL(max_th, h, k)
* avoid trouble by ignoring l = l_max
              do 70 l = l00, l_max, d_l
                if(l.eq.l_max) goto 70
                theta = ANGLE(h,k,l)
                call GET_F(f, S(h,k,l), l)
                tmp = INTENS(f, h, k, l, ok) * eps8
* find width of peak
                x = eps10
   80             if(.not.ok) goto 120
                  x = TWO * x
                  call GET_F(f, S(h,k,l+x), l+x)
                  if(INTENS(f, h, k, l+x, ok).gt.tmp .and.
     |                  x.le.eps2*d_l) goto 80
                if(.not.ok) goto 120
                l0 = max(l - x, ZERO)
                l1 = min(l + x, l_max)
                x = AGLQ16(h, k, l0, l1, ok)
                if(.not.ok) goto 110
*
* include weighting factors for radiation type
                if(rad_type.eq.X_RAY) then
                  x = TWO * x * W1(theta)
                else if(rad_type.eq.NEUTRN) then
                  x = TWO * x * W2(theta)
                else if(rad_type.eq.ELECTN) then
                  x = TWO * x * W3(theta)
                else
                  ok = .false.
                  goto 130
                endif
*
* see if not on l-axis
                if(.not.l_axis) then
* apply multiplicity factor
                  x = x * mltplcty
* if on boundary, apply appropriate weighting (mirror vs rotation only)
                  if(on_bndry) x = x * bnds_wt
                endif
                m = int(theta / d_theta) + 1
                if(finite_width) then
                  call CHWDTH(h,k,l0,l1,x,m,max_indx)
                else
                  spec(m) = spec(m) + x
                endif
   70         continue
            endif
          endif
   30   continue
   20 continue
      GETSPC = .true.
      return
  110 write(op,300) 'GLQ16 returned error in GETSPC.'
      return
  120 write(op,300) 'INTENS returned error in GETSPC.'
      return
  130 write(op,300) 'ERROR: Radiation type is undefined in GETSPC'
  999 return
  990 write(op,300) 'Illegal cell parameters in GETSPC.'
      write(op,250) '4*a0*b0-d0*d0 = ', FOUR * a0 * b0 - d0 * d0
      return
  200 format(1x, a, 2i4, 6x, 3a)
  250 format(1x, a, g12.5)
  300 format(1x, a)
  400 format(1x, i3)
      end
*
* ______________________________________________________________________
* Title: GET_BDS
* Author: MMJT
* Date: 1 July 1989
* Description:  This routine assigns reciprocal space vectors in the
* h,k plane within which integration can be confined. Weightings are
* assigned to off-axis spot intensities to allow for their
* multiplicity relative to spots on the 00l axis. Spots that occur on
* one of the boundaries are assigned special weighting depending on
* whether or not the boundary is also a mirror plane.
*
*      ARGUMENTS:
*           No arguments are used. All data is in 'COMMON'.
*
*      COMMON VARIABLES:
*            uses:      SymGrpNo, h_mirror, k_mirror
*
*        modifies:      h_start, k_start, h_end, k_end, mltplcty,
*                       bnds_wt, rot_only, pnt_grp
* ______________________________________________________________________
*
      subroutine GET_BDS()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
* 360 degrees, no symmetry (-1)
* note, the scan vectors are not used in this instance since there is
* an ambiguity between 0 and 360 degrees. We assign values anyway, so
* that the variables are at least initialized in a controlled way.
      if(SymGrpNo.eq.1) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ONE
        k_end    =  ZERO
        mltplcty =  ONE
        bnds_wt  =  ONE
        rot_only = .true.
      endif
*
* 180 degrees, rotation only (2/M, 1st setting)
      if(SymGrpNo.eq.2) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    = -ONE
        k_end    =  ZERO
        mltplcty =  TWO
        bnds_wt  =  ONE
        rot_only = .true.
      endif
*
* 180 degrees, vertical mirror (2/M, 2nd setting)
* we need to know which mirror plane.
      if(SymGrpNo.eq.3) then
        if(h_mirror .and. .not.k_mirror) then
          h_start  =  ONE
          k_start  =  ZERO
          h_end    = -ONE
          k_end    =  ZERO
          mltplcty =  TWO
          bnds_wt  =  HALF
          rot_only = .false.
        else if(k_mirror .and. .not.h_mirror) then
          h_start  =  ZERO
          k_start  =  ONE
          h_end    =  ZERO
          k_end    = -ONE
          mltplcty =  TWO
          bnds_wt  =  HALF
          rot_only = .false.
        else
* In theory, the following should never be needed, but just in case,
* let's bolster DIFFaX's immune system.
          write(op,400) 'DIFFaX is confused about vertical mirrors.'
          write(op,400) 'To be safe, symmetry is being set to -1'
          SymGrpNo = 1
          pnt_grp = '-1'
          h_start  =  ONE
          k_start  =  ZERO
          h_end    =  ONE
          k_end    =  ZERO
          mltplcty =  ONE
          bnds_wt  =  ONE
          rot_only = .true.
        endif
      endif
*
* 90 degrees, vertical mirrors (MMM)
      if(SymGrpNo.eq.4) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ZERO
        k_end    =  ONE
        mltplcty =  FOUR
        bnds_wt  =  HALF
        rot_only = .false.
      endif
*
* 120 degrees, rotation only (-3)
      if(SymGrpNo.eq.5) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    = -ONE
        k_end    =  ONE
        mltplcty =  THREE
        bnds_wt  =  ONE
        rot_only = .true.
      endif
*
* 60 degrees, vertical mirrors (-3M)
      if(SymGrpNo.eq.6) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ZERO
        k_end    =  ONE
        mltplcty =  SIX
        bnds_wt  =  HALF
        rot_only = .false.
      endif
*
* 90 degrees, rotation (4/M)
      if(SymGrpNo.eq.7) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ZERO
        k_end    =  ONE
        mltplcty =  FOUR
        bnds_wt  =  ONE
        rot_only = .true.
      endif
*
* 45 degrees, vertical mirrors (4/MMM)
      if(SymGrpNo.eq.8) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ONE
        k_end    =  ONE
        mltplcty =  EIGHT
        bnds_wt  =  HALF
        rot_only = .false.
      endif
*
* 60 degrees, rotation (6/M)
      if(SymGrpNo.eq.9) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ZERO
        k_end    =  ONE
        mltplcty =  SIX
        bnds_wt  =  ONE
        rot_only = .true.
      endif
*
* 30 degrees, vertical mirrors (6/MMM)
      if(SymGrpNo.eq.10) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ONE
        k_end    =  ONE
        mltplcty =  TWELVE
        bnds_wt  =  HALF
        rot_only = .false.
      endif
*
* integrate along 0 0 l only (axial)
* the following are somewhat arbitrary in this case. Assign values
* anyway just to make sure they are initialized
      if(SymGrpNo.eq.11) then
        h_start  =  ONE
        k_start  =  ZERO
        h_end    =  ONE
        k_end    =  ZERO
        mltplcty =  ONE
        bnds_wt  =  ONE
        rot_only = .true.
      endif
*
      return
  400 format(1x, a)
      end
*
* ______________________________________________________________________
* Title: GET_F
* Author: MWD and MMJT
* Date: 22 Mar 1989
* Description: This routine calculates the form factors for each layer.
* Since this routine is the main holdup for complex structures, it
* attempts to make use of shortcuts detected in OPTIMZ. The Debye-
* Waller diffuse background is assumed to be small and is not
* calculated in this version.
*
*      ARGUMENTS:
*            f   -  Array of layer form factors. (output).
*            Q2  -  Value of 1/d**2 at h,k,l. (input).
*            l   -  reciprocal lattice vector l-component. (input).
*
*      COMMON VARIABLES:
*            uses:  x_sf, n_sf, e_sf, rad_type, n_atoms, l_symmetry,
*                   one_B, l_n_atoms, a_type, hx_ky, a_pos, a_occup,
*                   a_B, l_actual, CENTRO, ELECTN, NEUTRN, X_RAY
*                   n_actual, n_layers
*
*        modifies:  no COMMON variables are modified
* ______________________________________________________________________
*
      subroutine GET_F(f, S2, l)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 S2, l
      complex*16 f(MAX_L)
*
      integer*4 i, j, m, n, type
      real*8 fact(MAX_TA), tmp(MAX_TA), tmp_sum, dot, e_factor, Q2
      parameter(e_factor = 0.023934D0)
      complex*16 ctmp(MAX_TA), f_uniq(MAX_L), ctmp_sum
*
      Q2 = QUARTER * S2
* Q2 = sin(theta)**2 / lamba**2
* determine scattering factors for each atom type
      if(rad_type.eq.X_RAY .or. rad_type.eq.ELECTN) then
        do 10 i = 1, n_atoms
* This empirical formula comes from p. 71 of
* "International Tables for X-ray Crystallography, Vol. IV"
* (The Kynoch Press: Birmingham, England), 1974.
          fact(i) = x_sf(1,i) * exp(-x_sf(2,i) * Q2) +
     |              x_sf(3,i) * exp(-x_sf(4,i) * Q2) +
     |              x_sf(5,i) * exp(-x_sf(6,i) * Q2) +
     |              x_sf(7,i) * exp(-x_sf(8,i) * Q2) +
     |              x_sf(9,i)
   10 continue
      else if(rad_type.eq.NEUTRN) then
        do 20 i = 1, n_atoms
          fact(i) = n_sf(i)
   20   continue
      endif
*
* get electron scattering factor from x-ray scattering factor
* s = 4 pi sin(theta) / lambda
* f_electron(s) = (8 pi**2 m e**2 / h**2) {Z - fx(s)} / s**2
*
*       = 0.023934 lambda**2 {Z - fx(s)} / sin(theta)**2
      if(rad_type.eq.ELECTN) then
        do 30 i = 1, n_atoms
          fact(i) = e_factor * ( dble(e_sf(i)) - fact(i) ) / Q2
   30   continue
      endif
*
      do 80 m = 1, n_actual
        tmp_sum  = ZERO
        ctmp_sum = C_ZERO
        do 35 n = 1, n_atoms
          tmp(n)  = ZERO
          ctmp(n) = C_ZERO
   35   continue
*
* First calculate the scattering factors of the unique layers.
* Check to see if f_uniq(m) will all be real and if Debye-Waller is
* invariant
* Note: hx_ky(j,m) contains h*a_pos(1,j,m) + k*a_pos(2,j,m)
*
        if(l_symmetry(m).eq.CENTRO .and. one_B(m)) then
          do 40 j = 1, l_n_atoms(m)
            type = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            tmp(type) = tmp(type) + a_occup(j,m) * cos(dot)
   40     continue
          do 45 j = 1, n_atoms
            tmp_sum = tmp_sum + tmp(j) * fact(j)
   45     continue
* NOTE: twice since centrosymmetric
          f_uniq(m) = TWO * exp(-a_B(1,m) * Q2) * tmp_sum
*
* Debye-Waller is not invariant
        else if(l_symmetry(m).eq.CENTRO) then
          do 50 j = 1, l_n_atoms(m)
            type = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            tmp(type) = tmp(type) + a_occup(j,m) *
     |              exp(-a_B(j,m) * Q2) * cos(dot)
   50     continue
          do 55 j = 1, n_atoms
            tmp_sum = tmp_sum + tmp(j) * fact(j)
   55     continue
* NOTE: twice since centrosymmetric
          f_uniq(m) = TWO * tmp_sum
*
* check if Debye-Waller is the only invariant
* f(i) will be complex
        else if(one_B(m)) then
          do 60 j = 1, l_n_atoms(m)
            type = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            ctmp(type) = ctmp(type) + a_occup(j,m) *
     |               dcmplx(cos(dot), sin(dot))
   60     continue
          do 65 j = 1, n_atoms
            ctmp_sum = ctmp_sum + ctmp(j) * fact(j)
   65     continue
          f_uniq(m) = exp(-a_B(1,m) * Q2) * ctmp_sum
*
* Nothing is invariant
        else
          do 70 j = 1, l_n_atoms(m)
            type = a_type(j,m)
            dot = hx_ky(j,m) + l*a_pos(3,j,m)
            ctmp(type) = ctmp(type) + a_occup(j,m) *
     |             exp(-a_B(j,m) * Q2) * dcmplx(cos(dot), sin(dot))
   70     continue
          do 75 j = 1, n_atoms
            ctmp_sum = ctmp_sum + ctmp(j) * fact(j)
   75     continue
          f_uniq(m) = ctmp_sum
        endif
   80 continue
*
* Now assign scattering factors to all the layers
      do 90 i = 1, n_layers
        f(i) = f_uniq(l_actual(i))
   90 continue
*
      return
      end
*
* ______________________________________________________________________
* Title: GET_G
* Author: MWD and MMJT
* Date: 18 Aug 1988; 15 Mar 1995
* Description: This function determines g_i, the a-priori probability
* that a layer of type i, will actually occur somewhere within the
* crystal.
* 'cnt' counts the l_alpha(i,i) = 1.0 terms. Then, l_g(i) = 1.0/cnt.
*
*      ARGUMENTS:
*            No input arguments.
*
*      COMMON VARIABLES:
*            uses:   n_layers, l_alpha
*
*        modifies:   l_g
*
*      GET_G returns logical .true. if all went well.
* ______________________________________________________________________
*
      logical function GET_G()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical singular, LUDCMP
      integer*4 i, j, cnt, index(MAX_L)
      real*8 sum, g_mat(MAX_L,MAX_L), Det
*
* external function
      external LUDCMP
* external subroutine (Some compilers need them declared external)
*      external LUBKSB
*
      GET_G = .false.
* set up matrix that represents the probabilities
* only n-1 equations are independent
      do 10 i = 1, n_layers - 1
        l_g(i) = ZERO
        sum = ZERO
        do 20 j = 1, n_layers
          sum = sum + l_alpha(j,i)
   20   continue
        sum = ONE / sum
* sum should actually be ONE
        do 30 j = 1, n_layers
          g_mat(i,j) = sum * l_alpha(i,j)
   30   continue
        g_mat(i,i) = g_mat(i,i) - ONE
   10 continue
      l_g(n_layers) = ONE
*
* the sum of the g's must be 1
      do 40 i = 1, n_layers
        g_mat(n_layers, i) = ONE
   40 continue
*
* before we invert the matrix, let's catch the pathological values
      cnt = 0
      do 50 i = 1, n_layers
        if(l_alpha(i,i).eq.ONE) cnt = cnt + 1
   50 continue
*
      if(cnt.ne.0) then
        singular = .true.
      else
        singular = .false.
      endif
*
      if(singular) then
        do 60 i = 1, n_layers
          if(l_alpha(i,i).eq.ONE) then
* arbitrarily assume such layers occur with equal probability
            l_g(i) = ONE / dble(cnt)
          else
            l_g(i) = ZERO
          endif
   60   continue
      else
* solve the matrix
        if(.not. LUDCMP(g_mat,index,n_layers,MAX_L,Det)) goto 100
        call LUBKSB(g_mat,l_g,index,n_layers,MAX_L)
      endif
*
* If we are here then all went well
      GET_G = .true.
*
* Are some of the layers non-existent?
      do 70 i = 1, n_layers
        if(l_g(i).lt.eps6) then
            write(op,410) 'WARNING: Layer ', i,
     |      ' does not occur in any significant quantity.'
          endif
   70 continue
*
      return
  100 write(op,400)
     | 'ERROR: Stacking probabilities give a singular matrix in GET_G'
      return
  400 format(1x, a)
  410 format(1x, a, i2, a)
      end
*
* ______________________________________________________________________
* Title: GET_MAT
* Author: MMJT
* Date: 21 Mar 1990; 14th July 1995; 21st July 1997
* Description:  This subroutine calculates the elements of 'mat', the
* stacking transition matrix ie. {alpha * exp(2*pi*u.R)}ij. The h and k
* components were pre-calculated in PRE_MAT. GET_MAT calculates the
* l-component to complete the matrix. However, the mat(i,i) terms
* have yet to have 1 subtracted from them. This is done in GET_S.
*
*      ARGUMENTS:
*            h   -  reciprocal vector h-component. (input).
*            k   -  reciprocal vector k-component. (input).
*            l   -  reciprocal lattice vector l-component. (input).
*
*      COMMON VARIABLES:
*            uses:  same_rz, n_layers, same_Bs, Bs_zero, mat1, c0, bc0,
*                   ca0, r_B33, r_B23, r_B31, l_r, PI2, there,
*                   fatsWalla_hk
*
*        modifies:  mat
* ______________________________________________________________________
*
      subroutine GET_MAT(h, k, l)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k
      real*8 l
*
      real*8 dot, twopi_l, fatsWaller
      integer*4 i, j
*
* set up matrix that represents the sequences
* Note: mat is in 'i,j' order.
      twopi_l = PI2 * l
      if(same_Bs) then
        if(all_Bs_zero) then
* all the fatsWalla_hk terms equal 1
          do 10 i = 1, n_layers
            do 20 j = 1, n_layers
              if(there(j,i)) then
                dot = twopi_l * l_r(3,j,i)
                mat(i,j) = mat1(i,j) * dcmplx( cos(dot), sin(dot) )
              else
                mat(i,j) = C_ZERO
              endif
   20       continue
   10     continue
        else
* all the fatsWalla terms are identical, but are less than 1
* fatsWalla_hk already contains the h-k component computed in PRE_MAT
          fatsWaller = fatsWalla_hk * exp(-(l*(QUARTER*a_B33*c0*l
     |         + HALF*(a_B23*bc0*k + a_B31*ca0*h))))
          do 30 i = 1, n_layers
            do 40 j = 1, n_layers
              if(there(j,i)) then
                dot = twopi_l * l_r(3,j,i)
                mat(i,j) = mat1(i,j) * fatsWaller
     |                   * dcmplx( cos(dot), sin(dot) )
              else
                mat(i,j) = C_ZERO
              endif
   40       continue
   30     continue
        endif
      else
* the fatsWalla terms differ. mat1 already contains the h-k component
        do 50 i = 1, n_layers
          do 60 j = 1, n_layers
            if(there(j,i)) then
              dot = twopi_l * l_r(3,j,i)
              if(Bs_zero(j,i)) then
                mat(i,j) = mat1(i,j) * dcmplx( cos(dot), sin(dot) )
              else
                mat(i,j) = mat1(i,j) * dcmplx( cos(dot), sin(dot) )
     |            * exp( -(l*(QUARTER*r_B33(j,i)*c0*l
     |            + HALF*(r_B23(j,i)*bc0*k + r_B31(j,i)*ca0*h))) )
              endif
            else
              mat(i,j) = C_ZERO
            endif
   60     continue
   50   continue
      endif
*
      return
      end
*
* ______________________________________________________________________
* Title: GET_S
* Author: MWD and MMJT
* Date: 5th Aug 1991
* Description:  This function determines the S's--the average scattered
* wave functions of each layer at h, k, l for crystals with an
* infinite number of layers. GET_MAT should be called prior to GET_S.
* Note, since GET_S is called billions of times from deep within the
* inner loops of DIFFaX's bowels, part of the matrix mat1 has been
* precalculated in PRE_MAT in order to speed things up a bit.
*
*      ARGUMENTS:
*            f   -  Array of layer form factors. (input).
*            s   -  Array of average layer wavefunctions. (output).
*            h   -  reciprocal vector h-component. (input).
*            k   -  reciprocal vector k-component. (input).
*            l   -  reciprocal lattice vector l-component. (input).
*
*      COMMON VARIABLES:
*            uses:  n_layers
*
*        modifies:  mat
*
*      GET_S returns logical .true. if all went well.
* ______________________________________________________________________
*
      logical function GET_S(f, s, h, k, l)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k
      real*8 l
      complex*16 f(MAX_L), s(MAX_L)
*
* i_ok is used by Linpack routines
      integer i_ok, index(MAX_L)
      integer*4 i
      complex*16 Det, s_tmp(2)
* external subroutines (Some compilers need them declared external)
* CGEFA and CGESL are Linpack routines
      external CGEFA, CGESL
*
      GET_S = .false.
*
* subtract identity matrix (need do this for diagonal terms only).
      do 10 i = 1, n_layers
        mat(i,i) = mat(i,i) - C_ONE
   10 continue
*
* Now solve the system of equations.
      if(n_layers.gt.2) then
* now call LINPACK routines
        call CGEFA(mat, MAX_L, n_layers, index, i_ok)
        if(i_ok.ne.0) goto 999
        call CGESL(mat, MAX_L, n_layers, index, s, 0)
      else if(n_layers.eq.2) then
* its a simple 2 x 2, so solve it directly
        Det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
        if(Det.eq.C_ZERO) goto 999
* copy s (remember, if infinitely thick, s = -f)
        s_tmp(1) = -s(1)
        s_tmp(2) = -s(2)
        s(1) = ( mat(1,2) * s_tmp(2) - mat(2,2) * s_tmp(1) ) / Det
        s(2) = ( mat(2,1) * s_tmp(1) - mat(1,1) * s_tmp(2) ) / Det
      else if(n_layers.eq.1) then
* only one layer, so solve it immediately
        s(1) = -f(1) / mat(1,1)
      endif
*
      GET_S = .true.
      return
  999 write(op,400) 'Solving for sequence produces a singular matrix.'
      write(op,401) h, k, l
      do 50 i = 1, n_layers
        write(op,402) i, f(i)
   50 continue
      return
  400 format(1x, 'GET_S:', a)
  401 format(1x, 'h = ', i3, ' k = ', i3, ' l = ', g12.5)
  402 format(5x, 'f(', i2, ') = (', g12.5, ',', g12.5, ')')
      end
*
* ______________________________________________________________________
* Title: GET_S2
* Author: MMJT
* Date: 5 Feb 1990
* Description:  This function determines the S's--the average scattered
* wave functions of each layer at h, k, l for crystals with only a
* finite number of layers, l_cnt. The equation being solved is
*
*   inv(Ident-T) * ((N+1)*Ident - inv(Ident-T)*(Ident-T**(N+1)) * F / N
*
*  where N = l_cnt, and T is the stacking probability matrix
*
*      ARGUMENTS:
*            f   -  Array of layer form factors. (input).
*            s   -  Array of average layer wavefunctions. (output).
*            h   -  reciprocal vector h-component. (input).
*            k   -  reciprocal vector k-component. (input).
*            l   -  reciprocal lattice vector l-component. (input).
*
*      COMMON VARIABLES:
*            uses:  n_layers, l_cnt
*
*        modifies:  mat
*
*      GET_S2 returns logical .true. if all went well.
* ______________________________________________________________________
*
      logical function GET_S2(f, s, h, k, l)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k
      real*8 l
      complex*16 f(MAX_L), s(MAX_L)
*
      logical ok, GET_S, MAT2N
      integer*4 i, j
      complex*16 ctmp, mat_n(MAX_L,MAX_L), tmp_mat(MAX_L,MAX_L)
* external functions
      external GET_S, MAT2N
*
      GET_S2 = .false.
*
* get matrix mat multiplied by itself l_cnt+1 times
      ok = MAT2N(mat_n)
      if(.not.ok) goto 990
*
* subtract identity matrix, and make a copy of mat.
      do 10 i = 1, n_layers
        do 20 j = 1, n_layers
          tmp_mat(j,i) = mat(j,i)
   20   continue
        mat_n(i,i) = mat_n(i,i) - C_ONE
   10 continue
*
* Multiply out -(Ident - T**(l_cnt+1))F.
      do 30 i = 1, n_layers
        ctmp = C_ZERO
        do 40 j = 1, n_layers
          ctmp = ctmp + mat_n(i,j) * f(j)
   40   continue
        s(i) = ctmp
   30 continue
* Next, solve. ie. inv(Ident - T) * (Ident - T**(l_cnt+1))*F
* where now s = (Ident - T**(l_cnt+1))*F
      ok = GET_S(f, s, h, k, l)
      if(.not.ok) goto 999
*
* Use result to build a new vector, and solve again.
* First, reconstruct mat.
      do 50 i = 1, n_layers
        s(i) = (s(i) - f(i) * dble(l_cnt + 1)) / dble(l_cnt)
        do 60 j = 1, n_layers
          mat(j,i) = tmp_mat(j,i)
   60   continue
   50 continue
* Solve with new RHS vector
      ok = GET_S(f, s, h, k, l)
      if(.not.ok) goto 999
*
      GET_S2 = .true.
      return
  990 write(op,400) 'MAT2N returned an error in GET_S2.'
      write(op,401) h, k, l
      return
  999 write(op,400) 'Solving for sequence produces a singular matrix.'
      write(op,401) h, k, l
      do 70 i = 1, n_layers
        write(op,402) i, f(i)
   70 continue
      return
  400 format(1x, 'GET_S2:', a)
  401 format(1x, 'h = ', i3, ' k = ', i3, ' l = ', g12.5)
  402 format(5x, 'f(', i2, ') = (', g12.5, ',', g12.5, ')')
      end
*
* ______________________________________________________________________
* Title: GET_SYM
* Author: MMJT
* Date: 19 June 1990
* Determines the symmetry of the diffraction pattern, thereby
* defining the smallest volume of reciprocal space which needs to
* be integrated over. There are only 10 kinematical diffraction point
* groups in the presence of streaking. Friedel's law ensures there
* is always a center of symmetry, and the possibility of streaking
* prohibits the cubic point groups. The point group symmetry number
* is returned as GET_SYM.
* The 10 groups are:
*
*              GET_SYM          point group
*           -------------------------------
*                 1:       |        -1
*                 2:       |        2/M(1) (diad along streaks)
*                 3:       |        2/M(2) (mirror contains streaks)
*                 4:       |        MMM
*                 5:       |        -3
*                 6:       |        -3M
*                 7:       |        4/M
*                 8:       |        4/MMM
*                 9:       |        6/M
*                10:       |        6/MMM
*
* The point group symbol is returned in the global character string
* 'pnt_grp'.
*
*      ARGUMENTS:
*            ok  -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:  cell_gamma, DoSymDump, cell_a, cell_b, PI, PI2
*                   RAD2DEG
*
*        modifies:  max_var, pnt_grp, h_mirror, k_mirror
*
*      GET_SYM returns one of the ten symmetry flags listed above.
* ______________________________________________________________________
*
      integer*4 function GET_SYM(ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical ok
*
      logical diad, triad, tetrad, hexad
      logical TST_ROT, TST_MIR
      logical cell90, cell120, eq_sides
      integer*4 idum, rot_sym
      real*8 tmp_var
*
* external functions
      external TST_ROT, TST_MIR
* external subroutine (Some compilers need them declared external)
*
* initialize random numbers in RAN3
      idum = -1
*
* initialize function
      GET_SYM = 0
*
      max_var = ZERO
      diad   = .false.
      triad  = .false.
      tetrad = .false.
      hexad  = .false.
      cell90 =  abs(cell_gamma - HALF*PI) .lt. HALF*PI*eps6
      cell120 = abs(cell_gamma - PI2/THREE) .lt. PI2*eps6/THREE
*
* sample reciprocal space to get an idea of the sort of intensities
* that are out there.
* test for vertical mirror (equivalent to a 2-fold, Friedel's law)
      tmp_var = max_var
      rot_sym = 2
      diad = TST_ROT(rot_sym, idum, ok)
      if(.not.ok) goto 997
      if(.not.diad) max_var = tmp_var
*
* if the cell angle is neither 90 nor 120 degrees, then the symmetry
* has to be either -1 or 2/M.
      if( .not.cell90 .and. .not.cell120 ) then
        if(DoSymDump) then
          write(sy,'(a)') ' '
          write(sy,220) cell_gamma * RAD2DEG
          write(sy,221)
        endif
*
        if(diad) then
          GET_SYM = 2
          pnt_grp = '2/M(1)'
        else
          GET_SYM = 1
          pnt_grp = '-1'
        endif
        goto 900
      endif
      eq_sides = abs(cell_a - cell_b) .le. HALF*eps6*(cell_a + cell_b)
      if(.not.eq_sides .and. DoSymDump) then
        write(sy,'(a)') ' '
        write(sy,225) cell_a, cell_b
        write(sy,226)
        write(sy,221)
      endif
*
* cell_gamma is either 90 or 120 degrees.
* if cell_a = cell_b, higher rotational symmetry is possible
      if(eq_sides) then
* Note, it is quite possible for an oblique cell (whose cell_gamma is
* not equal to 120 degrees) to have 3-fold symmetry. We do not test
* for this.
        tmp_var = max_var
        if(cell120) then
          rot_sym = 3
          triad = TST_ROT(rot_sym, idum, ok)
          if(.not.ok) goto 997
        else
          triad = .false.
        endif
        if(.not.triad) max_var = tmp_var
        hexad = diad.and.triad
        if(hexad.and.DoSymDump) write(sy,200)
        if(diad.and.cell90 .and. (.not.triad)) then
          tmp_var = max_var
          rot_sym = 4
          tetrad = TST_ROT(rot_sym, idum, ok)
          if(.not.ok) goto 997
          if(.not.tetrad) max_var = tmp_var
        else
          tetrad = .false.
        endif
      endif
*
* Now test for mirrors.
      tmp_var = max_var
      h_mirror = TST_MIR(1, idum, ok)
      if(.not.ok) goto 998
      if(.not.h_mirror) max_var = tmp_var
      tmp_var = max_var
      k_mirror = TST_MIR(2, idum, ok)
      if(.not.ok) goto 999
      if(.not.k_mirror) max_var = tmp_var
      tmp_var = max_var
      hk_mirror = TST_MIR(3, idum, ok)
      if(.not.ok) goto 999
      if(.not.hk_mirror) max_var = tmp_var
*
* If h_mirror does not equal k_mirror, then there cannot be a higher
* rotation symmetry than a diad. If, by some bizarre freak, this is
* inconsistent with the result of TST_ROT, choose the lower symmetry.
      if(h_mirror.neqv.k_mirror) then
        if(diad) then
          GET_SYM = 2
          pnt_grp = '2/M(1)'
        else
          GET_SYM = 3
          pnt_grp = '2/M(2)'
        endif
        goto 900
      endif
* Now check for combinations of mirrors and rotation axes.
*
* 6-fold
      if(hexad) then
        if(h_mirror.or.hk_mirror) then
          GET_SYM = 10
          pnt_grp = '6/MMM'
        else
          GET_SYM = 9
          pnt_grp = '6/M'
        endif
        goto 900
      endif
* 4-fold
      if(tetrad) then
        if(h_mirror.or.hk_mirror) then
          GET_SYM = 8
          pnt_grp = '4/MMM'
        else
          GET_SYM = 7
          pnt_grp = '4/M'
        endif
        goto 900
      endif
* 3-fold
      if(triad.and.(.not.diad)) then
        if(h_mirror.or.hk_mirror) then
          GET_SYM = 6
          pnt_grp = '-3M'
        else
          GET_SYM = 5
          pnt_grp = '-3'
        endif
        goto 900
      endif
* 2-fold
      if(diad) then
* handle special case of non-orthogonal mesh which has a diad,
* no triad, and one mirror. Diad prevails, vertical mirrors ignored.
        if((h_mirror.or.hk_mirror).and.cell90) then
          GET_SYM = 4
          pnt_grp = 'MMM'
        else
          GET_SYM = 2
          pnt_grp = '2/M(1)'
        endif
        goto 900
      endif
* if no symmetry has been detected opt for lowest symmetry
      GET_SYM = 1
      pnt_grp = '-1'
*
  900 return
  997 write(op,228) 'ERROR in GET_SYM: error returned by TST_ROT'
      write(op,229) '      while testing for ', rot_sym, '-fold axis.'
        return
  998 write(op,228) 'ERROR in GET_SYM: error returned by TST_MIR'
      write(op,228) '   while testing for mirror about the a - c plane'
        return
  999 write(op,228) 'ERROR in GET_SYM: error returned by TST_MIR'
      write(op,228) '   while testing for mirror about the b - c plane'
      return
  200 format(1x,'THE 2-FOLD AND 3-FOLD IMPLY 6-FOLD ROTATION SYMMETRY')
  220 format(1x, 'Cell_gamma = ', g12.5, ' degrees')
  221 format(1x, 'Rotational symmetry higher than 2-fold is unlikely.')
  225 format(1x, 'cell-a = ', g12.5, ' cell-b = ', g12.5)
  226 format(1x, 'Cell sides are not equal.')
  228 format(1x, a)
  229 format(1x, a, i2, a)
      end
*
* ______________________________________________________________________
* Title: GAUSSN
* Author: MMJT
* Date: 17 Feb 1990; 7 Mar 1995
* Description: This subroutine simulates Gaussian instrumental
* broadening. std_dev is in degrees. The algorithm used does not
* conserve intensity when std_dev is comparable to d_theta. Intensities
* at the extreme ends of the spectrum are corrupted slightly.
*
*      ARGUMENTS:
*            th2_low  -  lowest 2theta angle to consider. (input).
*
*      COMMON VARIABLES:
*            uses:  NONE, PI2, RAD2DEG, blurring, d_theta, FWHM
*                   th2_max
*
*        modifies:  brd_spc, spec
*
* ______________________________________________________________________
*
      subroutine GAUSSN(th2_low)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 th2_low
*
      integer*4 i, j, n_low, n_high, m
      real*8 k1, k2, k3, const, gss, std_dev, tmp, tmp1, tmp2
*
      if(FWHM.le.ZERO) goto 999
      std_dev = FWHM / sqrt(EIGHT * log(TWO))
*
* check that cut-off is reasonable
      if(th2_low.lt.ZERO .or. th2_low.ge.th2_max) then
        write(op,101) 'GAUSSN: Cut-off angle ', th2_low,
     |        ' is out of bounds. Angle reset to zero.'
        th2_low = ZERO
      endif
*
* th2_low is the angle relative to th2_min
* 2*d_theta is the angular step size
      n_low  = int(HALF*th2_low/d_theta) + 1
      n_high = int(HALF*(th2_max-th2_min)/d_theta) + 1
*
      const = TWO * RAD2DEG * d_theta
      k1 = const / (sqrt(PI2) * std_dev)
      k2 = HALF * (const / std_dev)**2
*
      do 10 i = 1, n_high
        brd_spc(i) = ZERO
   10 continue
*
* go out to 40 standard deviations, or to the end of the spectrum
      m = nint(TWO * TWENTY * std_dev / const)
      if(m.gt.n_high) m = n_high
      do 20 i = 0, m
        k3 = k2*dble(i*i)
        gss = k1*exp(-k3)
        do 30 j = n_low+1, n_high
          tmp1 = ZERO
          tmp2 = ZERO
          if((j-i).gt.n_low)  tmp1 = spec(j-i)
          if((j+i).le.n_high) tmp2 = spec(j+i)
          tmp = tmp1 + tmp2
          if(i.eq.0) tmp = HALF * tmp
          brd_spc(j) = brd_spc(j) + gss * tmp
   30   continue
   20 continue
      return
  999 write(op,101) 'Illegal FWHM ', FWHM, ' in GAUSSN.'
      write(op,100)'Gaussian instrumental broadening not added'
* kill blurring option
      blurring = NONE
      return
  100 format(1x, a)
  101 format(1x, a, g12.5, a)
      end
*
* ______________________________________________________________________
* Title: GLQ16
* Authors: MWD and MMJT
* Date: 6 April 1989; 13th July 95
*  This routine performs 16-point Gauss-Legendre quadrature on an
*  interval in reciprocal space. The interval is (h,k,a) to (h,k,b).
*  The routine calls INTENS at each of the 16 points, and if all goes
*  well, returns .true..
*  In the interests of speed, this routine has the option of calling
*  APPR_F, which returns interpolated f values for each of the 16
*  points. This modification slows down the procedure for structures
*  with few atoms per layer, but speeds it up if there are many atoms
*  per layer.
*
*      ARGUMENTS:
*            h   -  reciprocal lattice vector h-component. (input).
*            k   -  reciprocal lattice vector k-component. (input).
*            a   -  l-value of the lower bound of reciprocal
*                   lattice integration region. (input).
*            b   -  l-value of the upper bound of reciprocal
*                   lattice integration region. (input).
*            ok  -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:  a0, b0, c0, d0, recrsv
*
*        modifies:  intp_F
*
*      GLQ16 returns the integrated value.
* ______________________________________________________________________
*
      real*8 function GLQ16(h, k, a, b, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical ok
      integer*4 h, k
      real*8 a, b
*
      logical o, too_close
      real*8 INTENS, INTEN2
      real*8 c1, c2, x1, x2, x3, x4, x5, x6, x7, x8
      real*8         w1, w2, w3, w4, w5, w6, w7, w8
*
      parameter (x1 = 0.095012509837637440185D0)
      parameter (x2 = 0.281603550779258913230D0)
      parameter (x3 = 0.458016777657227386342D0)
      parameter (x4 = 0.617876244402643748447D0)
      parameter (x5 = 0.755404408355003033895D0)
      parameter (x6 = 0.865631202387831743880D0)
      parameter (x7 = 0.944575023073232576078D0)
      parameter (x8 = 0.989400934991649932596D0)
*
      parameter (w1 = 0.189450610455068496285D0)
      parameter (w2 = 0.182603415044923588867D0)
      parameter (w3 = 0.169156519395002538189D0)
      parameter (w4 = 0.149595988816576732081D0)
      parameter (w5 = 0.124628971255533872052D0)
      parameter (w6 = 0.095158511682492784810D0)
      parameter (w7 = 0.062253523938647892863D0)
      parameter (w8 = 0.027152459411754094852D0)
*
      integer*4 i, j
* f is approximated by a polynomial of order (n-1)
      integer*4 n
      parameter (n = 3)
      integer*4 list(n)
*
      real*8 Q2, l
      real*8 ag_l(16), samp_l(n)
      complex*16 f(MAX_L,16)
*
* external functions
      external INTENS, INTEN2
* external subroutines (Some compilers need them declared external)
*      external APPR_F, GET_F
*
* statement function
* Q2 is the value of 1/d**2 at hkl
      Q2(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
*
* initialize, to suppress compiler warnings
      GLQ16 = ZERO
*
* check that the integration range is legitimate
      if(b.lt.a) goto 999
* catch pathological values of a and b - i.e. zero integration range!
      if(b.eq.a) then
        ok = .true.
        goto 900
      endif
*
* let's interpolate ('hard-wired' in this version)
      intp_F = .true.
*
      c1 = HALF*(b - a)
      c2 = c1 + a
*
* set up the 16 l-values to sample
      ag_l(1) = -c1*x8 + c2
      ag_l(2) = -c1*x7 + c2
      ag_l(3) = -c1*x6 + c2
      ag_l(4) = -c1*x5 + c2
      ag_l(5) = -c1*x4 + c2
      ag_l(6) = -c1*x3 + c2
      ag_l(7) = -c1*x2 + c2
      ag_l(8) = -c1*x1 + c2
*
      ag_l( 9) = c1*x1 + c2
      ag_l(10) = c1*x2 + c2
      ag_l(11) = c1*x3 + c2
      ag_l(12) = c1*x4 + c2
      ag_l(13) = c1*x5 + c2
      ag_l(14) = c1*x6 + c2
      ag_l(15) = c1*x7 + c2
      ag_l(16) = c1*x8 + c2
*
      if(intp_F) then
*
* choose special values to sample (3 point interpolation)
        list(1) = 1
        list(2) = 8
        list(3) = 16
        samp_l(1) = ag_l(list(1))
        samp_l(2) = ag_l(list(2))
        samp_l(3) = ag_l(list(3))
*
* Deal with very rare cases when the spread in l values is too small
        too_close = (samp_l(1).eq.samp_l(3)) .or.
     |              (samp_l(1).eq.samp_l(2)) .or.
     |              (samp_l(2).eq.samp_l(3))
*
        if(.not.too_close) then
          call APPR_F(f, h, k, samp_l, ag_l, n, list, ok)
          if(.not.ok) goto 990
        else
* sample f values once, and set all 16 values over l-range to be equal
          call GET_F(f(1,1), Q2(h,k,ag_l(1)), ag_l(1))
          do 10 j = 1, n_layers
            do 20 i = 2, 16
              f(j,i) = f(j,1)
   20       continue
   10     continue
        endif
*
      else
* do not interpolate
        do 30 i = 1, 16
          call GET_F(f(1,i), Q2(h,k,ag_l(i)), ag_l(i))
   30   continue
*
      endif
*
* sum intensities
*
      o = .true.
      if(recrsv) then
*
      GLQ16 = c1 * (
     |w8*(INTENS(f(1,1),h,k,ag_l(1),o)+INTENS(f(1,16),h,k,ag_l(16),o))+
     |w7*(INTENS(f(1,2),h,k,ag_l(2),o)+INTENS(f(1,15),h,k,ag_l(15),o))+
     |w6*(INTENS(f(1,3),h,k,ag_l(3),o)+INTENS(f(1,14),h,k,ag_l(14),o))+
     |w5*(INTENS(f(1,4),h,k,ag_l(4),o)+INTENS(f(1,13),h,k,ag_l(13),o))+
     |w4*(INTENS(f(1,5),h,k,ag_l(5),o)+INTENS(f(1,12),h,k,ag_l(12),o))+
     |w3*(INTENS(f(1,6),h,k,ag_l(6),o)+INTENS(f(1,11),h,k,ag_l(11),o))+
     |w2*(INTENS(f(1,7),h,k,ag_l(7),o)+INTENS(f(1,10),h,k,ag_l(10),o))+
     |w1*(INTENS(f(1,8),h,k,ag_l(8),o)+INTENS(f(1, 9),h,k,ag_l( 9),o)))
*
      else
*
      GLQ16 = c1 * (
     |w8*(INTEN2(f(1,1),h,k,ag_l(1),o)+INTEN2(f(1,16),h,k,ag_l(16),o))+
     |w7*(INTEN2(f(1,2),h,k,ag_l(2),o)+INTEN2(f(1,15),h,k,ag_l(15),o))+
     |w6*(INTEN2(f(1,3),h,k,ag_l(3),o)+INTEN2(f(1,14),h,k,ag_l(14),o))+
     |w5*(INTEN2(f(1,4),h,k,ag_l(4),o)+INTEN2(f(1,13),h,k,ag_l(13),o))+
     |w4*(INTEN2(f(1,5),h,k,ag_l(5),o)+INTEN2(f(1,12),h,k,ag_l(12),o))+
     |w3*(INTEN2(f(1,6),h,k,ag_l(6),o)+INTEN2(f(1,11),h,k,ag_l(11),o))+
     |w2*(INTEN2(f(1,7),h,k,ag_l(7),o)+INTEN2(f(1,10),h,k,ag_l(10),o))+
     |w1*(INTEN2(f(1,8),h,k,ag_l(8),o)+INTEN2(f(1, 9),h,k,ag_l( 9),o)))
*
      endif
      ok = o
  900 return
  999 ok = .false.
      write(op,100) 'GLQ16: Illegal integration interval!'
      write(op,101) h, k, a, b
      return
  990 ok = .false.
      write(op,100) 'GLQ16: ERROR returned from APPR_F'
      write(op,101) h, k, a, b
      return
  100 format(1x, a)
  101 format(1x, 'h = ',i4,', k = ',i4,', l0 = ',g12.5,', l1 = ',g12.5)
      end
*
* ______________________________________________________________________
** Title: GOINTR
** Author: MMJT
** Date: 23 Oct 1989
** Description: This subroutine sets up the parameters necessary to
** integrate over an interval of reciprocal space.
**
**      ARGUMENTS:
**            ok  -  logical flag indicating all went well. (output).
**
**      COMMON VARIABLES:
**            uses:  cntrl, CFile
**
**        modifies:  no COMMON variables are modified
** ______________________________________________________________________
**
*      subroutine GOINTR(ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      logical ok
**
*      integer*4 i
*      real*8 AGLQ16, GLQ16
**
** external functions
*      external AGLQ16, GLQ16
** external subroutine (Some compilers need them declared external)
**      external INTEGR
**
*    1 write(op,100) ' '
*      write(op,100) 'INTEGRATING INTENSITY IN AN INTERVAL ALONG l. . .'
*      write(op,100) 'Enter 1 for adaptive quadrature.'
*      read(cntrl,*,err=1) i
*      if(CFile) write(op,101) i
*      if(i.eq.1) then
*        call INTEGR(AGLQ16,ok)
*      else
*        call INTEGR(GLQ16,ok)
*      endif
*      if(.not.ok) goto 999
*      i = 1
*    2 write(op,100) 'Enter 1 to integrate another interval.'
*      read(cntrl,*,err=2) i
*      if(CFile) write(op,101) i
*      if(i.eq.1) goto 1
**
*  999 return
*  100 format(1x, a)
*  101 format(1x, i3)
*      end
**
** ______________________________________________________________________
** Title: GOSADP
** Author: MMJT
** Date: 23 Oct 1989; 1st Mar 2000
** Description: This subroutine sets up a file whose name is given
** by 'outfile', which is derived from the input filename 'infile'.
** The selected area diffraction pattern (SADP) is calculated and is
** then written in binary format to 'outfile'.
**
**      ARGUMENTS:
**            infile   -  name of input file (input).
**            outfile  -  name of output file (output).
**            ok       -  logical flag indicating all went well.
**                                                      (output).
**
**      COMMON VARIABLES:
**            uses:  cntrl, CFile
**
**        modifies:  loglin, brightness
** ______________________________________________________________________
**
*      subroutine GOSADP(infile,outfile,ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      character*(*) infile, outfile
*      logical ok
**
*      integer*4 i_adapt, i_plane, io_err, hk_lim
*      real*8 AGLQ16, GLQ16, l_upper
**
** external functions
*      external AGLQ16, GLQ16
** external subroutines (Some compilers need them declared external)
**      external STREAK, GETFNM, GETSAD, WRTSAD
**
** Get a suitable file name, and open file.
*      call GETFNM(infile,outfile,'sadp',ok)
*      if(.not.ok) goto 999
** Open unformatted for binary write.
*      if(sa.ne.op) open(unit = sa, file = outfile, status = 'new',
*     |       form = 'unformatted', err = 990, iostat = io_err)
** some compilers need the following added for true binary output
**     |       ,recordtype = 'stream')
**
*      write(op,100) ' '
*      write(op,100)'CALCULATING SELECTED AREA DIFFRACTION PATTERN. . .'
**
* 1234 write(op,100) 'Enter 1 for adaptive quadrature.'
*      read(cntrl,*,err=1234) i_adapt
*      if(CFile) write(op,101) i_adapt
**
*    1 write(op,100) 'Choose a plane in reciprocal space to view.'
*      write(op,100) '       the l-axis is included by default.'
*      write(op,100) '1: k = 0.   2: h = 0.   3: h = k.   4: h = -k.'
*      read(cntrl,*,err=1) i_plane
*      if(CFile) write(op,101) i_plane
*      if(i_plane.lt.1 .or. i_plane.gt.4) then
*        if(CFile) then
*          write(op,100) 'Illegal choice. Default is k = 0.'
*          i_plane = 1
*        else
*          write(op,100) 'Illegal choice. Choose again.'
*          goto 1
*        endif
*      endif
**
** Get upper bounds. The SADP will be square, so we only need one bound.
** Choose l-axis, since this is common to all 4 views and makes scaling
** simpler if more than one SAD pattern is requested.
* 2345 write(op,100) 'Enter maximum value of l.'
*      read(cntrl,*,err=2345) l_upper
*      if(CFile) write(op,104) l_upper
**
** 8-bit images or 16-bit images?
* 3456 write(op,100) 'Choose the bit-depth of the image.'
*      write(op,100) '   8: - 8-bits  16: - 16-bits'
*      read(cntrl,*,err=3456) bitdepth
*      if(CFile) write(op,101) bitdepth
*      if(bitdepth.ne.8 .and. bitdepth.ne.16) then
*        write(op,100) 'Illegal bit-depth.'
*        if(CFile) then
*          write(op,100) 'Using 8-bit as the default'
*          bitdepth = 8
*        else
*          write(op,100) 'Re-enter. . .'
*          goto 3456
*        endif
*      endif
*      if(bitdepth.eq.16) then
** Bypass issue of signed or unsigned format. Use range 0 - 32767 only.
*        write(op,100) 
*     | 'File will be saved in unsigned 16-bit (big-endian) format.'
*        maxsad = FIFTEENBITS - ONE
*      else
*        write(op,100) 'File will be saved in unsigned 8-bit format.'
*        maxsad = EIGHTBITS - ONE
*      endif
**
** Logarithmic or linear intensity scaling?
*    2 write(op,100) 'Choose intensity scaling'
*      write(op,100) '   0: - Logarithmic  1: - Linear'
*      read(cntrl,*,err=2) loglin
*      if(CFile) write(op,101) loglin
*      if(loglin.lt.0 .or. loglin.gt.1) then
*        write(op,100) 'Illegal intensity scaling type.'
*        if(CFile) then
*          write(op,100) 'Using linear as the default'
*          loglin = 1
*        else
*          write(op,100) 'Re-enter. . .'
*          goto 2
*        endif
*      endif
**
** By how much do we wish to saturate?
*    3 write(op,100) 'Enter a brightness (must be +ve)'
*      if(bitdepth.eq.16) then
*        write(op,100)'  1 - scales intensities to the range 0 - 32767'
*        write(op,100)' 10 - scales intensities to the range 0 - 327670,'
*        write(op,100)'      but values above 65535 will be saturated'
*      else
*        write(op,100)'  1 - scales intensities to the range 0 - 255'
*        write(op,100)' 10 - scales intensities to the range 0 - 2550,'
*        write(op,100)'      but values above 255 will be saturated'
*      endif
*      read(cntrl,*,err=3) brightness
*      if(CFile) write(op,104) brightness
*      if(brightness.le.ZERO) then
*        write(op,100) 'Illegal value for brightness. Must be positive'
*        if(CFile) then
*          if(loglin.eq.0) then
*            write(op,100) 'Using default of 1 for logarithmic scaling'
*            brightness = ONE
*          else
*            write(op,100) 'Using default of 10 for linear scaling'
*            brightness = TEN
*          endif
*        else
*          if(loglin.eq.0) then
*            write(op,100) '1 is a good default for logarithmic scaling'
*          else
*            write(op,100) '10 is a good default for linear scaling'
*          endif
*          write(op,100) 'Re-enter. . .'
*          goto 3
*        endif
*      endif
**
** Generate intensity data for the SAD pattern.
*      if(i_adapt.eq.1) then
*        call GETSAD(AGLQ16, i_plane, l_upper, hk_lim, infile, ok)
*      else
*        call GETSAD( GLQ16, i_plane, l_upper, hk_lim, infile, ok)
*      endif
**
*      if(ok) call WRTSAD(outfile, i_plane, l_upper, hk_lim, ok)
**
*  999 return
*  990 write(op,102) 'Problems opening output file ', outfile
*      write(op,103) 'IOSTAT = ', io_err
*      ok = .false.
*      return
*  100 format(1x, a)
*  101 format(1x, i3)
*  102 format(1x, 2a)
*  103 format(1x, a, i5)
*  104 format(1x, g12.5)
*      end
**
** ______________________________________________________________________
** Title: GOSPEC
** Author: MMJT
** Date: 23 Oct 1989
** Description: This subroutine sets up a file whose name is given
** by 'outfile', which is derived from the input filename 'infile'.
** The powder pattern data is then written to 'outfile', which is closed
** by the subroutine 'WRTSPC'.
**
**      ARGUMENTS:
**            infile   -  name of input file (input)
**            outfile  -  name of output file (output)
**            ok       -  logical flag indicating all went well.
**                                                      (output).
**
**      COMMON VARIABLES:
**            uses:  blurring, CFile, GAUSS, LORENZ, NONE, PS_VGT
**                   PV_GSS, PV_LRN, cntrl
**
**        modifies:  full_brd
** ______________________________________________________________________
**
*      subroutine GOSPEC(infile,outfile,ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      character*(*) infile, outfile
*      logical ok
**
*      logical GETSPC, TRMSPC, RDRNGE
*      integer*4 io_err
*      real*8 AGLQ16, GLQ16, cut_off
**
** external functions
*      external RDRNGE, GETSPC, AGLQ16, GLQ16, TRMSPC
** external subroutines (Some compilers need them declared external)
**      external GETFNM, PV, GAUSSN, LORNZN, WRTSPC
**
*      write(op,100) ' '
*      write(op,100) 'CALCULATING POWDER DIFFRACTION PATTERN. . .'
**
** get angular range and step
*      ok = RDRNGE()
*      if(.not.ok) goto 999
** create a new filename to save spectrum data to
*      call GETFNM(infile, outfile, 'spc', ok)
*      if(.not.ok) goto 999
*      if(sp.ne.op) open(unit = sp, file = outfile, status = 'new',
*     |            err = 990, iostat = io_err)
*      full_brd = 1
* 3456 write(op,100) 'Enter 1 for adaptive quadrature on broad peaks.'
*      read(cntrl,*,err=3456) full_brd
*      if(CFile) write(op,101) full_brd
*      if(full_brd.eq.1) then
*        ok = GETSPC(AGLQ16, infile)
*      else
*        ok = GETSPC(GLQ16, infile)
*      endif
** suppress the huge peak near the origin if required
*      cut_off = ZERO
*      if(ok.and.(th2_min.eq.ZERO).and.trim_origin) ok = TRMSPC(cut_off)
*      if(ok) then
*        if(blurring.eq.GAUSS) then
*          write(op,104) 'Gaussian'
*          call GAUSSN(cut_off)
*        else if(blurring.eq.LORENZ) then
*          write(op,104) 'Lorentzian'
*          call LORNZN(cut_off)
*        else if(blurring.eq.PS_VGT) then
*          write(op,104) 'pseudo-Voigt'
*          call PV(cut_off)
*        else if(blurring.eq.PV_GSS) then
*          write(op,104) 'Gaussian'
*          call PV(cut_off)
*        else if(blurring.eq.PV_LRN) then
*          write(op,104) 'Lorentzian'
*          call PV(cut_off)
*        else if(blurring.ne.NONE) then
*          write(op,100)
*     |         'Instrumental broadening type is undefined in GOSPEC.'
*        endif
*      endif
*      if(ok) call WRTSPC(outfile, ok)
**
*  999 return
*  990 write(op,102) 'Problems opening output file ', outfile
*      write(op,103) 'IOSTAT = ', io_err
*      ok = .false.
*      return
*  100 format(1x, a)
*  101 format(1x, i3)
*  102 format(1x, 2a)
*  103 format(1x, a, i5)
*  104 format(1x, 'Adding ', a, ' instrumental broadening . . .')
*      end
**
** ______________________________________________________________________
** Title: GOSTRK
** Author: MMJT
** Date: 23 Oct 1989
** Description: This subroutine sets up a file whose name is given
** by 'outfile', which is derived from the input filename 'infile'.
** The streak intensity is then written to 'outfile', which is closed
** by the subroutine 'STREAK'.
**
**      ARGUMENTS:
**            infile   -  name of input file. (input).
**            outfile  -  name of output file. (output).
**            ok       -  logical flag indicating all went well.
**                                                       (output).
**
**      COMMON VARIABLES:
**            uses:  cntrl, CFile
**
**        modifies:  no COMMON variables are modified
** ______________________________________________________________________
**
*      subroutine GOSTRK(infile,outfile,ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      character*(*) infile, outfile
*      logical ok
**
*      integer*4 i, io_err
*      real*8 AGLQ16, GLQ16
**
** external functions
*      external AGLQ16, GLQ16
** external subroutines (Some compilers need them declared external)
**      external STREAK, GETFNM
**
*      call GETFNM(infile, outfile, 'str', ok)
*      if(.not.ok) goto 999
*      if(sk.ne.op) open(unit = sk, file = outfile, status = 'new',
*     |                              err = 990, iostat = io_err)
* 4567 write(op,100) ' '
*      write(op,100) 'CALCULATING INTENSITY ALONG A STREAK. . .'
*      write(op,100) 'Enter 1 for adaptive quadrature'
*      read(cntrl,*,err=4567) i
*      if(CFile) write(op,101) i
** select integration function
*      if(i.eq.1) then
*        call STREAK(AGLQ16, outfile, ok)
*      else
*        call STREAK( GLQ16, outfile, ok)
*      endif
**
*  999 return
*  990 write(op,102) 'Problems opening output file ', outfile
*      write(op,103) 'IOSTAT = ', io_err
*      ok = .false.
*      return
*  100 format(1x, a)
*  101 format(1x, i3)
*  102 format(1x, 2a)
*  103 format(1x, a, i5)
*      end
**
** ______________________________________________________________________
* Title: HKL_LIM
* Author: MMJT
* Date: 2 August 1989
* Obtains upper limits of h, k, and l given the global variable
* 'max_angle' (in radians).
* The limits are returned in the global variables h_bnd, k_bnd and
* l_bnd. HKL_LIM may need to decrease the value of lambda if h_bnd
* and k_bnd are too small to allow adequate off-axis symmetry testing.
* lambda is restored in OPTIMZ after symmetry testing.
*
*      ARGUMENTS:
*           No arguments are used. All data is in 'COMMON'.
*
*      COMMON VARIABLES:
*            uses:  a0, b0, c0, d0, lambda
*
*        modifies:  h_bnd, k_bnd, l_bnd, lambda
* ______________________________________________________________________
*
      subroutine HKL_LIM()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 HKLUP
      real*8 y
*
* HKLUP returns the maximum value of h, k or l given 'max_angle'
      HKLUP(y) =  int(TWO * sin(HALF*max_angle) / (lambda*sqrt(y)))
*
* define upper h, k, l values consistent with max_angle
    1 h_bnd = HKLUP(a0)
      k_bnd = HKLUP(b0)
      l_bnd = HKLUP(c0)
*
* make sure bounds are not too small. This could occur for a small
* unit cell, a small value of max_angle, or a long wavelength.
      if(h_bnd.lt.2 .or. k_bnd.lt.2) then
        lambda = HALF * lambda
        goto 1
      endif
*
* Make sure bounds are not too large either
      if(h_bnd.gt.10) h_bnd = 10
      if(k_bnd.gt.10) k_bnd = 10
      if(l_bnd.gt.10) l_bnd = 10
*
      return
      end
*
* ______________________________________________________________________
* Title: INTEGR
* Author: MMJT and MWD
* Date: 15 Feb 1990; 7 Mar 1995; 28th May 1996
* Description:  This routine integrates intensity from
*               h,k,l0 to h,k,l1.
*
*      ARGUMENTS:
*            FN   -  Function name passed by reference. The
*                    choice is between GLQ16 (non-adaptive
*                    Gauss-Legendre integration), and AGLQ16
*                    (adaptive Gauss-Legendre integration). (input).
*            ok   -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:   a0, b0, c0, d0, lambda, cntrl, CFile, xplcit,
*                    X_RAY, rad_type, th2_max
*
*        modifies:   no COMMON variables are modified
* ______________________________________________________________________
*
      subroutine INTEGR(FN, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 FN
      logical ok
*
      logical divided
      integer*4 h, k
      real*8 l0, l1, max_th, x, W4, ANGLE, theta, l, S, Q2
      real*8 t1, sum, tmp, LL, fact, d_th, l_tmp
*
* external function, passed by reference
      external FN
* external subroutines (Some compilers need them declared external)
*      external XYPHSE, PRE_MAT, GET_F
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
* ANGLE is the Bragg angle (in radians) of the h,k,l plane
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
* LL is the maximum l value for a given h,k
      LL(theta,h,k) = sqrt((Q2 * sin(theta) * sin(theta)
     |                    - h*h*a0 - k*k*b0 - h*k*d0)/ c0)
* W4 is the X-ray polarization factor
      W4(theta) = HALF * (ONE + (cos(TWO*theta))**2)
*
      max_th = HALF * th2_max
      fact = TWO / lambda
      Q2 = fact * fact
* set increment to a safe value, in case l1-l0 is too large
      d_th = eps3 * DEG2RAD
*
   10 write(op,400) 'Enter h, k, l0, l1'
      read(cntrl,*,err=10)  h, k, l0, l1
      if(CFile) write(op,401) h, k, l0, l1
* check values
      if(l1.eq.l0) then
        write(op,400) 'Illegal input: l1 equals l0'
        if(CFile) then
          l1 = l0 + d_th
          write(op,403) 'l1 is set to ', l1
        else
          goto 10
        endif
      endif
* make sure we are not going to blow up at the origin
      if(h.eq.0 .and. k.eq.0 .and. rad_type.eq.ELECTN) then
        if(l0*l1.le.ZERO) then
          write(op,400)
     |   'Cannot integrate across the origin for electron radiation'
          write(op,400) 'Re-enter. . .'
          goto 10
        endif
      endif
* Finally, check angles are meaningful
      if(S(h,k,l0).gt.Q2 .or. S(h,k,l1).gt.Q2) then
        if(S(h,k,l0).gt.Q2) write(op,402) h, k, l0,
     |            ' exceeds 180 degree scattering angle!'
        if(S(h,k,l1).gt.Q2) write(op,402) h, k, l1,
     |            ' exceeds 180 degree scattering angle!'
        goto 10
      endif
* get angles corresponding to start and stop
* check if we need to break the integration region into two parts
* because h,k,l, and h,k,-l may subtend the same angle
      divided  = .false.
      if(l0.le.ZERO .and. l1.le.ZERO) then
* use Friedel's law, and keep l +ve. Swap l0 and l1
        h = -h
        k = -k
        tmp = -l0
        l0 = -l1
        l1 = tmp
      else if(l0.lt.ZERO .and. l1.gt.ZERO) then
        h = -h
        k = -k
        l_tmp = l1
        l1 = -l0
        l0 = ZERO
        divided = .true.
      endif
* swap if in reverse order
      if(l0.gt.l1) then
        tmp = l0
        l0 = l1
        l1 = tmp
      endif
*
      sum = ZERO
   30 max_th = ANGLE(h,k,l1)
      t1     = ANGLE(h,k,l0)
      l1 = l0
      call XYPHSE(h, k)
      call PRE_MAT(h, k)
* integrate each d_th's range of reciprocal space
      do 40 theta = t1, max_th-eps14, d_th
        l0 = l1
        tmp = min(d_th, max_th-theta)
        l1 = LL(theta+tmp,h,k)
        x = FN(h,k,l0,l1,ok)
        if(.not.ok) goto 999
        if(rad_type.eq.X_RAY) x = x * W4(theta + HALF*tmp)
        sum = sum + x
   40 continue
* do we need to integrate the other part?
      if(divided) then
* goto -h,-k,-l and continue
        h = -h
        k = -k
        l0 = ZERO
        l1 = l_tmp
        divided = .false.
        goto 30
      endif
      write(op,403) 'Integrated intensity = ', sum
  999 return
  400 format(1x, a)
  401 format(1x, 2i3, 2g12.5)
  402 format(1x, 2i3, g12.5, a)
  403 format(1x, a, g15.8)
      end
*
* ______________________________________________________________________
* Title: INTEN2
* Author: MMJT
* Date: 8 Apr 1992
* Description: This routine determines the intensity at (h,k,l)
* in reciprocal space, for an explicitly defined stacking sequence.
*
*      ARGUMENTS:
*            f   -  Layer form factors. (input).
*            h   -  reciprocal lattice vector h-component. (input).
*            k   -  reciprocal lattice vector k-component. (input).
*            l   -  reciprocal lattice vector l-component. (input).
*            ok  -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:  same_layer, l_seq, l_r, n_layers, l_phi, l_cnt
*                   PI2
*
*        modifies:  wavefn
*
*      INTEN2 returns the intensity at h, k, l
* ______________________________________________________________________
*
      real*8 function INTEN2(f, h, k, l, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical ok
      integer*4 h, k
      real*8 l
      complex*16 f(MAX_L)
*
      integer*4 i, j, m
      real*8 twopi_l, dot, tmp
      complex*16 phi(MAX_L, MAX_L), z, z_to_n
*
      twopi_l = PI2 * l
*
* 'ok' is not used, but is included for compatibility with INTENS
      ok = .true.
*
* Is there only one layer? If so, let's get it over with.
      if(l_cnt.eq.1) then
        wavefn = f(l_seq(1))
        goto 900
      endif
*
* Check for an obvious optimization when all layers are identical.
      if(same_layer) then
        i = l_seq(1)
        dot = PI2*(h*l_r(1,i,i) + k*l_r(2,i,i) + l*l_r(3,i,i))
        z = dcmplx( cos(dot), sin(dot) )
        tmp = dot / PI2
* check we are not about to execute 0.0 / 0.0
        if(abs(tmp - nint(tmp)) .le. eps5) then
          wavefn = f(i) * l_cnt
        else
* sum the series
          dot = dot * l_cnt
          z_to_n = dcmplx( cos(dot), sin(dot) )
          wavefn = f(i) * (C_ONE - z_to_n) / (C_ONE - z)
        endif
        goto 900
      endif
*
* Else do it the long way
* Get phases
      do 10 i = 1, n_layers
        do 20 j = 1, n_layers
          dot = twopi_l * l_r(3,j,i)
          phi(j,i) = l_phi(j,i) * dcmplx( cos(dot), sin(dot) )
   20   continue
   10 continue
* Count down to the first layer (we know l_cnt is greater than 1)
* Initialize wavefunction to the scattering factor of the last layer
      wavefn = f(l_seq(l_cnt))
      do 30 m = l_cnt - 1, 1, -1
        i = l_seq(m)
        j = l_seq(m+1)
        wavefn = f(i) + wavefn * phi(j,i)
   30 continue
*
* Normalize to the number of layers
  900 INTEN2 = wavefn * conjg(wavefn) / l_cnt
*
      return
      end
*
* ______________________________________________________________________
* Title: INTENS
* Author: MWD and MMJT
* Date: 10 Feb 1989
* Description: This routine determines the intensity at (h,k,l)
* in reciprocal space, for recursive stacking. For this function
* to be called, 'rcrsv' must be TRUE.
* Note: The diffuse background is handled automatically by the
* recursion algorithm.
*
*      ARGUMENTS:
*            f   -  Layer form factors. (input).
*            h   -  reciprocal lattice vector h-component. (input).
*            k   -  reciprocal lattice vector k-component. (input).
*            l   -  reciprocal lattice vector l-component. (input).
*            ok  -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:  only_real, l_g, n_layers, inf_thick
*
*        modifies:  no COMMON variables are modified
*
*      INTENS returns the intensity at h, k, l
* ______________________________________________________________________
*
      real*8 function INTENS(f, h, k, l, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical ok
      integer*4 h, k
      real*8 l
      complex*16 f(MAX_L)
*
      logical GET_S, GET_S2
      integer*4 i
      real*8 sum, x
      complex*16 s(MAX_L)
*
* external function
      external GET_S, GET_S2
* external subroutines (Some compilers need them declared external)
*      external GET_MAT
*
      sum = ZERO
      call GET_MAT(h, k, l)
      if(inf_thick) then
* initialize s to -f, since mat is -(1 - T)
        do 10 i = 1, n_layers
          s(i) = - f(i)
   10   continue
        ok = GET_S(f, s, h, k, l)
      else
* s is initialized inside GET_S2, from where GET_S is called
        ok = GET_S2(f, s, h, k, l)
      endif
      if(ok) then
* only use real part of f(i) if that's all that's there
        if(only_real) then
          do 20 i = 1, n_layers
            sum = sum + l_g(i) * dble(f(i)) * dble(s(i))
   20     continue
          sum = TWO * sum
          do 30 i = 1, n_layers
            x = dble(f(i))
            sum = sum - l_g(i) * x*x
   30     continue
        else
* must use complex part of f(i)
          do 40 i = 1, n_layers
            sum = sum + l_g(i) * dble(conjg(f(i)) * s(i))
   40     continue
          sum = TWO * sum
          do 50 i = 1, n_layers
            x = abs(f(i))
            sum = sum - l_g(i) * x*x
   50     continue
        endif
      endif
*
      INTENS = sum
*
      return
      end
*
* ______________________________________________________________________
** Title: LAYCNT
** Author: MMJT
** Date: 3 Oct 1989
** Description: Counts the number of integer arguments in the character
** string 'line', and returns the answer in LAYCNT. Legal separators
** are, blanks, tabs, and commas. If LAYCNT encounters anything other
** than a digit (0 - 9), or one of the three legal separators, LAYCNT
** is set to -1. This routine is similar to CNTARG, accept it will only
** accept integers as arguments.
**
**      ARGUMENTS:
**            line  -  Line of characters to be parsed. (input).
**
**      LAYCNT returns the number of integer arguments in the line.
**      Returns -1 if an illegal character was detected.
** ______________________________________________________________________
**
*      integer*4 function LAYCNT(line)
*      implicit none
**
*      character*(*) line
**
*      logical in_num, legit
*      integer*4 blank, tab, comma, nought, nine, i, j, num_cnt, lin_len
*      parameter (tab=9, blank=32, comma=44, nought=48, nine=57)
**
*      in_num = .false.
*      legit = .true.
*      num_cnt = 0
*      LAYCNT = -1
**
*      lin_len = len(line)
*      do 10 i = 1, lin_len
*        j = ichar(line(i:i))
*        if(j.eq.tab .or. j.eq.blank .or. j.eq.comma) then
*          in_num = .false.
*        else if(j.ge.nought .and. j.le.nine) then
*          if(.not.in_num) then
*            in_num = .true.
*            num_cnt = num_cnt + 1
*          endif
*        else
** illegal character
*          legit = .false.
*        endif
*   10 continue
**
*      if(legit) LAYCNT = num_cnt
**
*      return
*      end
*
* ______________________________________________________________________
* Title: LENGTH
* Author: MMJT
* Date: 20 Nov 1989
* Description:  This function returns the length of the first group
* of contiguous, non-space characters in the passed string. If the
* string has no blanks, the declared length of the string is returned.
*
*      ARGUMENTS:
*            string  -  Character string whose length is needed.
*                                                            (input).
*      LENGTH returns the string length.
* ______________________________________________________________________
*
      integer*4 function LENGTH(string)
      implicit none
*
      character string*(*)
*
      integer*4 i
*
      i = index(string,' ')
      if(i.eq.0) then
        LENGTH = len(string)
      else
        LENGTH = i-1
      endif

      return
      end
*
* ______________________________________________________________________
* Title: LORNZN
* Author: MMJT
* Date: 17 Feb 1990; 7 Mar 1995
* Description: This subroutine performs the Lorentzian
* instrumental broadening. FWHM is in degrees. Does not conserve
* intensity well when FWHM is comparable to d_theta. Data at the
* extreme ends of the spectrum are corrupted slightly.
*
*      ARGUMENTS:
*            th2_low  -  lowest 2theta angle to consider. (input).
*
*      COMMON VARIABLES:
*            uses:  th2_max, d_theta, FWHM, spec, NONE, PI, RAD2DEG
*                   blurring
*
*        modifies:  brd_spc
* ______________________________________________________________________
*
      subroutine LORNZN(th2_low)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 th2_low
*
      integer*4 i, j, n_low, n_high
      real*8 k1, k2, k3, const, lrnz, tmp, tmp1, tmp2
*
      if(FWHM.le.ZERO) goto 999
* check that cut-off is reasonable
      if(th2_low.lt.ZERO .or. th2_low.ge.th2_max) then
        write(op,101) 'LORNZN: Cut-off angle ', th2_low,
     |        ' is out of bounds. Angle reset to zero.'
        th2_low = ZERO
      endif
*
* th2_low is the angle relative to th2_min
* 2*d_theta is the angular step size
      n_low  = int(HALF*th2_low/d_theta) + 1
      n_high = int(HALF*(th2_max-th2_min)/d_theta) + 1
*
      const = TWO * RAD2DEG * d_theta
      k1 = const * TWO / (PI * FWHM)
      k2 = (const * TWO / FWHM)**2
*
      do 10 i = 1, n_high
        brd_spc(i) = ZERO
   10 continue
*
      do 20 i = 0, n_high
        k3 = ONE + k2*i*i
        lrnz = k1 / k3
        do 30 j = n_low+1, n_high
          tmp1 = ZERO
          tmp2 = ZERO
          if((j-i).gt.n_low)  tmp1 = spec(j-i)
          if((j+i).le.n_high) tmp2 = spec(j+i)
          tmp = tmp1 + tmp2
          if(i.eq.0) tmp = HALF * tmp
          brd_spc(j) = brd_spc(j) + lrnz * tmp
   30   continue
   20 continue
      return
  999 write(op,101) 'Illegal FWHM ', FWHM, ' in LORNZN()'
      write(op,100) 'Lorentzian instrumental broadening not added'
* kill blurring option
      blurring = NONE
      return
  100 format(1x, a)
  101 format(1x, a, g12.5, a)
      end
*
* ______________________________________________________________________
* Title: LUBKSB
* Author: MWD, adapted from
* "Numerical Recipes: The Art of Scientific Computing."
* Date: 18 Aug 1988
*  Description: Solves the set of linear equations ax = b,
*  where a, x and b contain real variables.
*  Here, a is input as the LU-decomposition of a, determined by the
*  routine LUDCMP. index is input as the permutation vector returned
*  by LUDCMP. b is input as the right hand side vector, and returns
*  with the solution x. a, n, MAX_N and index are not modified by this
*  routine. In DIFFaX, LUDCMP and LUBKSB are used to solve for l_g(i),
*  the a-priori probability that layer i exists within the crystal.
*
*      ARGUMENTS:
*            a      -  LU-decomposed square matrix of real numbers.
*                                                           (input).
*            b      -  vector of real numbers, the right hand side of
*                      ax = b, is input. The solution x is output.
*            index  -  vector containing the record of the row
*                      permutation effected by the partial pivoting
*                      carried out in CLUDCM (input).
*            n      -  size of the square matrix. (input).
*            MAX_N  -  physical dimension of a (MAX_N x MAX_N).(input).
* ______________________________________________________________________
*
      subroutine LUBKSB(a,b,index,n,MAX_N)
      include 'DIFFaX.par'
*     save
*
      integer*4 n, MAX_N, index(MAX_N)
      real*8 a(MAX_N,MAX_N), b(MAX_N)
*
      integer*4 i, i2, j, row
      real*8 sum
*
      i2 = 0
      do 20 i = 1, n
        row = index(i)
        sum = b(row)
        b(row) = b(i)
        if(i2.ne.0) then
          do 10 j = i2, i-1
            sum = sum - a(i,j) * b(j)
   10     continue
        else if(abs(sum).ne.ZERO) then
          i2 = i
        endif
        b(i) = sum
   20 continue
      do 40 i = n, 1, -1
        sum = b(i)
        do 30 j = i+1, n
          sum = sum - a(i,j) * b(j)
   30   continue
        b(i) = sum / a(i,i)
   40 continue
      return
      end
*
* ______________________________________________________________________
* Title: LUDCMP
* Author: MWD, adapted from
* "Numerical Recipes: The Art of Scientific Computing."
* Date: 18 Aug 1988
*  Description: This is an LU decomposition routine, and accepts
*  real*8 variables.
*  Given an n x n matrix a, with physical dimension MAX_N, this
*  routine replaces it by the LU decomposition of a rowwise permutation
*  of itself. a and n are input. a is the LU decomposed output; index
*  is an output vector which records the row permutation affected by
*  the partial pivoting; Det is the determinant of a. This routine is
*  used in combination with LUBKSB to solve linear equations. In DIFFaX,
*  these routines are used to solve for l_g(i), the a-priori
*  probability that layer i exists within the crystal.
*  LUDCMP returns .false. if the matrix turns out to be singular.
*
*      ARGUMENTS:
*            a      -  Square matrix of real numbers to LU-decompose is
*                      input. a is then replaced by the
*                      LU-decomposed result.
*            index  -  output vector which records the row permutation
*                      affected by the partial pivoting. (output).
*            n      -  size of the square matrix. (input).
*            MAX_N  -  physical dimension of a (MAX_N x MAX_N).(input).
*            Det    -  determinant of the matrix. (output).
*
*      LUDCMP returns logical .true. if all proceeded happily.
* ______________________________________________________________________
*
      logical function LUDCMP(a,index,n,MAX_N,Det)
      include 'DIFFaX.par'
*     save
*
      integer*4 n, MAX_N, index(MAX_N)
      real*8 a(MAX_N,MAX_N), Det
*
      integer*4 L_MAX, i, j, m, row
      parameter (L_MAX = 100)
      real*8 tiny, tmp(L_MAX), sum, max, tmp2
      parameter (tiny = 1.0D-20)
*
      LUDCMP = .false.
      Det = ONE
      if(n.gt.L_MAX) then
        write(op,400) 'Matrix too large for LUDCMP'
        return
      endif
      do 20 i = 1, n
          max = ZERO
        do 10 j = 1, n
          if(abs(a(i,j)).gt.max) max = abs(a(i,j))
   10   continue
        if(max.eq.ZERO) goto 200
        tmp(i) = ONE / max
   20 continue
      do 90 j = 1, n
        do 40 i = 1, j-1
          sum = a(i,j)
          do 30 m = 1, i-1
            sum = sum - a(i,m) * a(m,j)
   30     continue
          a(i,j) = sum
   40   continue
        max = ZERO
        do 60 i = j, n
          sum = a(i,j)
          do 50 m = 1, j-1
            sum = sum - a(i,m) * a(m,j)
   50     continue
          a(i,j) = sum
          tmp2 = tmp(i)*abs(sum)
          if(abs(tmp2).ge.max) then
            row = i
            max = tmp2
          endif
   60   continue
        if(j.ne.row) then
          do 70 m = 1, n
            tmp2 = a(row,m)
            a(row,m) = a(j,m)
            a(j,m) = tmp2
   70     continue
          Det = -Det
          tmp(row) = tmp(j)
        endif
        index(j) = row
        if(abs(a(j,j)).eq.ZERO) a(j,j) = tiny
        tmp2 = ONE / a(j,j)
        do 80 i = j+1, n
          a(i,j) = a(i,j) * tmp2
   80   continue
        Det = Det * a(j,j)
   90 continue
      LUDCMP = .true.
      return
  200 continue
      return
  400 format(1x, a)
      end
*
* ______________________________________________________________________
* Title: L_STEP
* Author: MMJT
* Date: 13 Aug 1989
* Description:  L_STEP attempts to determine whether or not there
* are any sharp peaks, and if so, their spacing. The algorithm inspects
* the sum of certain permutations of the stacking vector z_components.
* The components are chosen so that their values should be independent
* of the layer origins chosen by the user. ie. Rz(i,i) and
* Rz(i,j) + Rz(j,i) and Rz(i,j) + Rz(j,k) + Rz(k,i), etc...
* In the interests of clarity, the permutations of the indices are
* written out explicity, instead of being generated by a subroutine
* (which would considerably shorten the code). We stop searching after
* combinations of 4 stacking vectors, since then the structure is too
* complicated to be attempting to trick DIFFaX into believing that we
* have found sharp peaks. Under these circumstances, work out the
* diffraction pattern the long, but reliable, way.
*
*      ARGUMENTS:
*            ok  -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:  n_layers, there, l_r
*
*        modifies:  any_sharp
*
*      L_STEP returns the l-increment that sharp peaks are likely to
*      be found at.
* ______________________________________________________________________
*
      real*8 function L_STEP(ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical ok
*
      real*8 tmp, z_step
      logical YRDSTK, resonant, decided
      integer*4 i1, i2, i3, i4
*
* external function
      external YRDSTK
*
* initialize return value
      L_STEP = ZERO
*
      resonant = .true.
      decided  = .false.
      z_step   = ZERO
*
* Check z-components of Rii stacking vectors
* if any of the transitions do not exist, set resonant to false.
      do 10 i1 = 1, n_layers
        if(resonant) then
          if(there(i1,i1)) then
            decided = .true.
            tmp = l_r(3,i1,i1)
            resonant = resonant .and. YRDSTK(z_step, tmp, ok)
            if(.not.ok) goto 990
          endif
        endif
   10 continue
*
* Rii terms do not occur (ie. no layer i will stack to another layer i)
* We must therefore check z-components of Rij + Rji sequences (i.ne.j).
      if((n_layers.gt.1) .and. .not.decided) then
        do 20 i1 = 1, n_layers
          do 30 i2 = i1 + 1, n_layers
            if(resonant) then
              if(there(i2,i1).and.there(i1,i2)) then
                decided = .true.
                tmp = l_r(3,i2,i1) + l_r(3,i1,i2)
                resonant = resonant.and.YRDSTK(z_step, tmp, ok)
                if(.not.ok) goto 991
              endif
            endif
   30   continue
   20 continue
      endif
*
* No Rij + Rji sequences occur.
* Check z-components of Rij + Rjk + Rki sequences (where i.ne.j.ne.k).
      if((n_layers.gt.2) .and. .not.decided) then
        do 40 i1 = 1, n_layers
          do 50 i2 = i1 + 1, n_layers
            do 60 i3 = i2 + 1, n_layers
              if(resonant) then
* There are 2 permutations
                if(there(i2,i1).and.
     |             there(i3,i2).and.
     |             there(i1,i3)) then
                   decided = .true.
                   tmp = l_r(3,i2,i1) + l_r(3,i3,i2) + l_r(3,i1,i3)
                   resonant = resonant .and. YRDSTK(z_step, tmp, ok)
                   if(.not.ok) goto 992
                endif
                if(there(i3,i1).and.
     |             there(i2,i3).and.
     |             there(i1,i2).and.resonant) then
                   decided = .true.
                   tmp = l_r(3,i3,i1) + l_r(3,i2,i3) + l_r(3,i1,i2)
                   resonant = resonant .and. YRDSTK(z_step, tmp, ok)
                   if(.not.ok) goto 993
                endif
              endif
   60       continue
   50     continue
   40   continue
      endif
*
* No Rij + Rjk + Rki sequences occur.
* Check z-components of Rij + Rjk + Rkl + Rli sequences
* (where i.ne.j.ne.k.ne.l).
      if((n_layers.gt.3) .and. .not.decided) then
        do 70 i1 = 1, n_layers
          do 80 i2 = i1 + 1, n_layers
            do 90 i3 = i2 + 1, n_layers
              do 100 i4 = i3 + 1, n_layers
                if(resonant) then
* There are 6 permutations
                  if(there(i2,i1).and.
     |               there(i3,i2).and.
     |               there(i4,i3).and.
     |               there(i1,i4)) then
                       decided = .true.
                       tmp = l_r(3,i2,i1) + l_r(3,i3,i2) +
     |                       l_r(3,i4,i3) + l_r(3,i1,i4)
                       resonant = resonant.and.YRDSTK(z_step,tmp,ok)
                       if(.not.ok) goto 994
                  endif
                  if(there(i2,i1).and.
     |               there(i4,i2).and.
     |               there(i3,i4).and.
     |               there(i1,i3).and.resonant) then
                       decided = .true.
                       tmp = l_r(3,i2,i1) + l_r(3,i4,i2) +
     |                        l_r(3,i3,i4) + l_r(3,i1,i3)
                       resonant = resonant.and.YRDSTK(z_step,tmp,ok)
                       if(.not.ok) goto 995
                  endif
                  if(there(i3,i1).and.
     |               there(i2,i3).and.
     |               there(i4,i2).and.
     |               there(i1,i4).and.resonant) then
                       decided = .true.
                       tmp = l_r(3,i3,i1) + l_r(3,i2,i3) +
     |                       l_r(3,i4,i2) + l_r(3,i1,i4)
                       resonant = resonant.and.YRDSTK(z_step,tmp,ok)
                       if(.not.ok) goto 996
                  endif
                  if(there(i3,i1).and.
     |               there(i4,i3).and.
     |               there(i2,i4).and.
     |               there(i1,i2).and.resonant) then
                       decided = .true.
                       tmp = l_r(3,i3,i1) + l_r(3,i4,i3) +
     |                       l_r(3,i2,i4) + l_r(3,i1,i2)
                       resonant = resonant.and.YRDSTK(z_step,tmp,ok)
                       if(.not.ok) goto 997
                  endif
                  if(there(i4,i1).and.
     |               there(i2,i4).and.
     |               there(i3,i2).and.
     |               there(i1,i3).and.resonant) then
                       decided = .true.
                       tmp = l_r(3,i4,i1) + l_r(3,i2,i4) +
     |                       l_r(3,i3,i2) + l_r(3,i1,i3)
                       resonant = resonant.and.YRDSTK(z_step,tmp,ok)
                       if(.not.ok) goto 998
                  endif
                  if(there(i4,i1).and.
     |               there(i3,i4).and.
     |               there(i2,i3).and.
     |               there(i1,i2).and.resonant) then
                       decided = .true.
                       tmp = l_r(3,i4,i1) + l_r(3,i3,i4) +
     |                       l_r(3,i2,i3) + l_r(3,i1,i2)
                       resonant = resonant.and.YRDSTK(z_step,tmp,ok)
                       if(.not.ok) goto 999
                  endif
                endif
  100         continue
   90       continue
   80     continue
   70   continue
      endif
*
* If there is no stacking sequence that can bring us back to a layer
* similar to that at the origin after 4 layers, then this
* structure is sufficiently complicated that we may be better
* off doing adaptive integration anyway. (d_l still equals 0.0.)
*
      if(decided.and.resonant .and. (tmp.ne.ZERO)) then
        L_STEP = ONE / tmp
        any_sharp = .true.
      else
        L_STEP = ZERO
        any_sharp = .false.
      endif
*
      return
  990 write(op,200)
      write(op,201) i1, i1
      return
  991 write(op,202)
      write(op,203) i1, i2, i2, i1
      return
  992 write(op,202)
      write(op,204) i1, i2, i2, i3, i3, i1
      return
  993 write(op,202)
      write(op,204) i1, i3, i3, i2, i2, i1
      return
  994 write(op,202)
      write(op,205) i1, i2, i2, i3, i3, i4, i4, i1
      return
  995 write(op,202)
      write(op,205) i1, i2, i2, i4, i4, i3, i3, i1
      return
  996 write(op,202)
      write(op,205) i1, i3, i3, i2, i2, i4, i4, i1
      return
  997 write(op,202)
      write(op,205) i1, i3, i3, i4, i4, i2, i2, i1
      return
  998 write(op,202)
      write(op,205) i1, i4, i4, i2, i2, i3, i3, i1
      return
  999 write(op,202)
      write(op,205) i1, i4, i4, i3, i3, i2, i2, i1
      return
  200 format(1x,'L_STEP: Non-physical z-component of stacking vector.')
  201 format(1x,'Rz(',i2,',',i2,') = 0.0')
  202 format(1x,'L_STEP:Non-physical z-components of stacking vectors')
  203 format(1x,'Rz(',i2,',',i2,') + Rz(',i2,',',i2,') = 0.0')
  204 format(1x,'Rz(',i2,',',i2,')',2(' + Rz(',i2,',',i2,')'),' = 0.0')
  205 format(1x,'Rz(',i2,',',i2,')',3(' + Rz(',i2,',',i2,')'),' = 0.0')
      end
*
* ______________________________________________________________________
* Title: MATMUL
* Author: MMJT
* Date: 5 Feb 1990
* Description:  This subroutine multiplies the complex matrices
* 'a' and 'b', of logical size n x n. Result is returned in 'a'.
*
*      ARGUMENTS:
*            a   -  Complex*16 array to store result. (input and output).
*            b   -  Complex*16 array. (input).
*            n   -  Logical size of matrices. (input).
*
* ______________________________________________________________________
*
      subroutine MATMUL(a, b, n)
      include 'DIFFaX.par'
*     save
*
      integer*4 n
      complex*16 a(MAX_L,MAX_L), b(MAX_L,MAX_L)
*
      integer*4 i, j, m
      complex*16 c(MAX_L,MAX_L), ctmp
*
* first copy a into c
      do 10 j = 1, n
        do 20 i = 1, n
          c(i,j) = a(i,j)
   20   continue
   10 continue
*
      do 30 j = 1, n
        do 40 i = 1, n
          ctmp = C_ZERO
          do 50 m = 1, n
            ctmp = ctmp + c(i,m) * b(m,j)
   50     continue
          a(i,j) = ctmp
   40   continue
   30 continue
*
      return
      end
*
* ______________________________________________________________________
* Title: MAT2N
* Author: MMJT
* Date: 5 Feb 1990
* Description:  This function multiplies the matrix 'mat' by itself
* n = l_cnt+1 times. In order to speed up the process, n has been broken
* down into its binary representation (by the function BINPOW, which is
* called once in OPTIMZ), and the result is given by
*
*         mat**n = mat**n0 * mat**n1 * mat**n2 * mat**n3 * etc
*
* where ni = 2**i
* and n = n0 + n1 + n2 + n3 + . . .
*
* mat**ni is given by (mat**n(i-1)) * (mat**n(i-1))
* ie. mat**8 = (mat**4) * (mat**4)
* and similarly mat**4 = (mat**2) * (mat**2)
*
* n must be such that n <= RCSV_MAX+1 <= 2**(MAX_BIN+1) - 1
*
*      ARGUMENTS:
*            a   -  Complex*16 array to store result. (output).
*
*      COMMON VARIABLES:
*            uses:  n_layers, pow, max_pow
*
*        modifies:  No COMMON variable are modified.
*
*      MAT2N returns TRUE if all went well.
* ______________________________________________________________________
*
      logical function MAT2N(a)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      complex*16 a(MAX_L,MAX_L)
*
      integer*4 i, j
      complex*16 tmp_mat(MAX_L,MAX_L,MAX_BIN)
*
* external subroutine (Some compilers need them declared external)
      external MATSQR, MATMUL
*
      MAT2N = .false.
*
* copy mat into the first 2-dimensional tmp_mat array. Initialize a
* to be the identity matrix.
      do 20 j = 1, n_layers
        do 30 i = 1, n_layers
          a(i,j) = C_ZERO
          tmp_mat(i,j,1) = mat(i,j)
   30   continue
        a(j,j) = C_ONE
   20 continue
*
      do 40 i = 1, max_pow-1
        if(pow(i).eq.1) then
          call MATMUL(a, tmp_mat(1,1,i), n_layers)
        endif
        call MATSQR(tmp_mat(1,1,i+1), tmp_mat(1,1,i), n_layers)
   40 continue
      if(pow(max_pow).eq.1)
     |         call MATMUL(a, tmp_mat(1,1,max_pow), n_layers)
*
      MAT2N = .true.
*
      return
      end
*
* ______________________________________________________________________
* Title: MATSQR
* Author: MMJT
* Date: 5 Feb 1990
* Description:  This subroutine multiplies the complex matrix
* 'b', of logical size n x n, by itself. Result is returned in 'a'.
*
*      ARGUMENTS:
*            a   -  Complex*16 array to store result. (output).
*            b   -  Complex*16 array to be 'squared'. (input).
*            n   -  Logical size of matrices. (input).
*
* ______________________________________________________________________
*
      subroutine MATSQR(a, b, n)
      include 'DIFFaX.par'
*     save
*
      integer*4 n
      complex*16 a(MAX_L,MAX_L), b(MAX_L,MAX_L)
*
      integer*4 i, j, m
      complex*16 ctmp
*
      do 10 j = 1, n
        do 20 i = 1, n
          ctmp = C_ZERO
          do 30 m = 1, n
            ctmp = ctmp + b(i,m) * b(m,j)
   30     continue
          a(i,j) = ctmp
   20   continue
   10 continue
*
      return
      end
*
* ______________________________________________________________________
** Title: NAMER
** Author: MMJT
** Date: 13 Oct 1988
** Description: This subroutine creates the appropriate filenames for
** the various files needed.
** It generates the filename by reading the structure data file name.
** It scans the name to see if a period ('.') is already within
** the name (ie. if the data file is called 'Structure.dat'). It deletes
** the characters after the period and appends whatever is in append.
** 'name1' and 'name2' are assumed to be the same physical length.
** This subroutine is called by GETFNM
**
**      ARGUMENTS:
**            name1   -  The name of the input data file. (input).
**            name2   -  A derivative filename. (output).
**            append  -  A token to be appended to name2. (input).
** ______________________________________________________________________
**
*      subroutine NAMER(name1,name2,append)
*      include 'DIFFaX.par'
**     save
**
*      character*(*) name1, name2, append
**
*      integer*4 applen, namlen, i, idot
*      integer*4 LENGTH
**
** external function
*      external LENGTH
**
** get length of the string holding the filename
*      namlen = len(name1)
*      name2  = ' '
*      applen = LENGTH(append)
**
*      idot = index(name1,'.')
*      if(idot.eq.0) idot = index(name1,' ')
*      if(idot.eq.0) then
*        if(namlen.lt.MAX_NAM) then
*          idot = namlen + 1
*        else
*          idot = MAX_NAM
*          namlen = MAX_NAM
*        endif
*      endif
** truncate root filename so that the appendage will appear correctly
*      if(idot+applen.gt.namlen) idot = namlen - applen
**
*      do 10 i = 1, idot-1
*        write(name2(i:i),'(a)') name1(i:i)
*   10 continue
**
*      write(name2(idot:idot),'(a)') '.'
**
*      do 20 i = 1, applen
*        write(name2((idot+i):(idot+i)),'(a)') append(i:i)
*   20 continue
**
*      do 30 i = idot+applen+1, namlen
*        write(name2(i:i),'(a)') ' '
*   30 continue
**
*      return
*      end
**
** ______________________________________________________________________
* Title: NMCOOR
* Author: MWD
* Date: 18 Aug 1988
* Description:  This subroutine multiplies the relative coordinates
* of each atom by 2*pi, which is the useful form for calculating phases.
*
*      ARGUMENTS:
*            No input arguments.
*
*      COMMON VARIABLES:
*            uses:  n_actual, l_n_atoms, a_pos, PI2
*
*        modifies:  a_pos
* ______________________________________________________________________
*
      subroutine NMCOOR()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 i, j
*
      do 10 i = 1, n_actual
        do 20 j = 1, l_n_atoms(i)
          a_pos(1,j,i) = a_pos(1,j,i) * PI2
          a_pos(2,j,i) = a_pos(2,j,i) * PI2
          a_pos(3,j,i) = a_pos(3,j,i) * PI2
   20   continue
   10 continue
*
      return
      end
*
* ______________________________________________________________________
** Title: NXTARG
** Author: MMJT
** Date: 21 JUL 1997
** Description: Removes the next argument from the character string
** 'line', and returns that argument in the character string 'arg'.
**
**      ARGUMENTS:
**            line  -  Line of characters to be parsed. (input).
**            arg   -  First argument found in 'line'.  (output).
**
**      Returns  1 if an argument was found.
**      Returns  0 if no argument was found.
**      Returns <0 if an error occurred.
** ______________________________________________________________________
**
*      integer*4 function NXTARG(line,arg)
*      implicit none
*      save
**
*      character*(*) line, arg
**
*      character*200 tmpline
*      integer*4 blank, tab, comma, i, j, i1, i2, lin_len
*      parameter (tab=9, blank=32, comma=44)
**
*      NXTARG = 0
**
*      lin_len = len(line)
**
*      i = 0
** Strip leading blanks
*   10 i = i + 1
*        if(i.gt.lin_len) goto 997
*        j = ichar(line(i:i))
*        if(j.eq.tab .or. j.eq.blank .or. j.eq.comma) goto 10
*      i1 = i
**
** Read off the non-blank characters
*   20 i = i + 1
*        if(i.gt.lin_len) goto 30
*        j = ichar(line(i:i))
*        if(j.eq.tab .or. j.eq.blank .or. j.eq.comma) then
*          NXTARG = 1
*          goto 30
*        else
*          goto 20
*        endif
*   30 continue
**
*      i2 = i
*      if(i2.lt.i1) goto 998
*      write(arg, '(a)', err=999) line(i1:i2)
*      write(tmpline, '(a)', err=999) line(i2:lin_len)
*      write(line, '(a)', err=999) tmpline(1:lin_len)
**
*      return
*  997 NXTARG = -3
*      return
*  998 NXTARG = -2
*      return
*  999 NXTARG = -1
*      return
*      end
**
** ______________________________________________________________________
* Title: OPTIMZ
* Author: MWD and MMJT
* Date: 8 Apr 1992; 15 Mar 1995; 24 July 1997
* Description:  This routine determines if any shortcuts can be taken
* in the calculations.
*
*      ARGUMENTS:
*            rootnam  -  The name of the input data file. (input).
*            ok       -  logical flag indicating all went well.
*                                                      (output).
*
*      COMMON VARIABLES:
*            uses:  a0, b0, d0, n_layers, l_alpha, l_actual, l_n_atoms,
*                   a_B, r_B11, r_B22, r_B33, r_B12, r_B23, r_B31,
*                   l_symmetry, l_seq, DoSymDump, SymGrpNo, lambda,
*                   l_cnt, CENTRO, PI, l_r, n_actual, finite_width
*
*        modifies:  there, one_B, Bs_zero, same_Bs, only_real, l_rz,
*                   same_layer, max_var, no_trials, tolerance, SymGrpNo,
*                   lambda, check_sym, has_l_mirror, theta1, theta2
*                   h_end, h_start, k_end, k_start, pnt_grp, same_rz
*                   xplcit, formfactor, all_Bs_zero,
*                   a_B11,a_B22,a_B33,a_B12,a_B23,a_B31
* ______________________________________________________________________
*
      subroutine OPTIMZ(rootnam, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      character*(*) rootnam
      logical ok
*
      character*31 sym_fnam
      logical EQUALB, BINPOW, GET_G
      integer*4 GET_SYM, i, j, j2, m, n, LENGTH
      real*8 HKANGL, h_val, k_val
      real*8 x, error, tmp, incr, z, old_lambda
      logical did_it(MAX_L,MAX_L)
*
* external functions
      external GET_SYM, LENGTH, EQUALB, BINPOW, GET_G
* external subroutines (Some compilers need them declared external)
*      external GETFNM, GET_BDS, CHK_SYM, NMCOOR, OVERLP
*
* statement function
* HKANGL is the angle between the vector (h_val,k_val,0) and (1,0,0)
      HKANGL(k_val,h_val) = atan2(k_val*sqrt(a0*b0 - d0*d0*QUARTER),
     |                              (h_val*a0 + k_val*d0*HALF))
*
* set up logic table for stacking transitions
      do 10 i = 1, n_layers
        do 20 j = 1, n_layers
          there(j,i) = l_alpha(j,i).ge.eps7
   20   continue
   10 continue
*
* see if there are any overlapping atoms
      call OVERLP()
*
* multiply all atom coordinates by 2*PI
      call NMCOOR()
*
* If calculation is to be recursive for a finite number of layers,
* then store the binary form of l_cnt+1 in an array for efficient
* matrix multiplication.
      ok = GET_G()
      if(recrsv .and. .not.inf_thick) then
        ok = BINPOW(l_cnt+1)
        if(.not.ok) then
          write(op,202) 'ERROR returned by BINPOW to OPTIMZ'
          write(op,201) 'The argument passed was l_cnt+1 = ', l_cnt+1
          goto 999
        endif
      endif
*
* see if Debye-Waller coefficients are same for all atoms in a layer
      do 30 i = 1, n_layers
        x = ZERO
        j2 = l_actual(i)
        m = l_n_atoms(j2)
        do 40 j = 1, m
          x = x + a_B(j,j2)
   40   continue
        x = x / dble(m)
        error = ZERO
* find absolute deviation of coefficients
        do 50 j = 1, m
          error = error + abs(a_B(j,j2) - x)
   50   continue
* get relative error
        if(x.ne.ZERO) error = error / (x * m)
        one_B(j2) = abs(error).le.eps3
   30 continue
*
* Check that the layer uncertainty factors are physically reasonable.
      do 60 i = 1, n_layers
        do 65 j = 1, n_layers
          if(there(j,i)) then
* check on r_B12
            x = r_B11(j,i)*r_B22(j,i)*a0*b0 -
     |          r_B12(j,i)*r_B12(j,i)*ab0*ab0
            if(x.lt.ZERO) then
              write(op,500) 'C12'
              write(op,501) i, j
              write(op,502) 'C11 and C22'
              ok = .false.
              goto 999
            endif
* check on r_B23
            x = r_B22(j,i)*r_B33(j,i)*b0*c0 -
     |          r_B23(j,i)*r_B23(j,i)*bc0*bc0
            if(x.lt.ZERO) then
              write(op,500) 'C23'
              write(op,501) i, j
              write(op,502) 'C22 and C33'
              ok = .false.
              goto 999
            endif
* check on r_B31
            x = r_B11(j,i)*r_B33(j,i)*a0*c0 -
     |          r_B31(j,i)*r_B31(j,i)*ca0*ca0
            if(x.lt.ZERO) then
              write(op,500) 'C13'
              write(op,501) i, j
              write(op,502) 'C11 and C33'
              ok = .false.
              goto 999
            endif
          endif
   65   continue
   60 continue
*
* see if stacking 'uncertainty' coefficients are same for all layers.
* Special flag if they are all zero.
      all_Bs_zero = .true.
      do 70 i = 1, n_layers
        do 80 j = 1, n_layers
          Bs_zero(j,i) =
     |         r_B11(j,i).eq.ZERO .and. r_B22(j,i).eq.ZERO .and.
     |         r_B33(j,i).eq.ZERO .and. r_B12(j,i).eq.ZERO .and.
     |         r_B23(j,i).eq.ZERO .and. r_B31(j,i).eq.ZERO
          all_Bs_zero = all_Bs_zero .and. Bs_zero(j,i)
   80   continue
   70 continue
*
* run through all 6 coefficients
* until it is clear that some are different
                  same_Bs = EQUALB(r_B11, a_B11)
      if(same_Bs) same_Bs = EQUALB(r_B22, a_B22)
      if(same_Bs) same_Bs = EQUALB(r_B33, a_B33)
      if(same_Bs) same_Bs = EQUALB(r_B12, a_B12)
      if(same_Bs) same_Bs = EQUALB(r_B23, a_B23)
      if(same_Bs) same_Bs = EQUALB(r_B31, a_B31)
*
* see if all layers are centrosymmetric
      only_real = .true.
      do 90 i = 1, n_actual
        only_real = only_real .and. (l_symmetry(i).eq.CENTRO)
   90 continue
*
* see if all z-components of the stacking vectors are the same
      l_rz = ZERO
      n = 0
      do 100 i = 1, n_layers
        do 110 j = 1, n_layers
          if(there(j,i)) then
            l_rz = l_rz + l_r(3,j,i)
            n = n + 1
          endif
  110   continue
  100 continue
      l_rz = l_rz / dble(n)
      error = ZERO
      do 120 i = 1, n_layers
        do 130 j = 1, n_layers
          if(there(j,i)) error = error + abs(l_r(3,j,i) - l_rz)
  130   continue
  120 continue
      same_rz = abs(error).le.eps4
*
* If the stacking is explicit, check to see if all the layers are
* the same
      same_layer = .false.
      if(xplcit) then
        if(l_cnt.eq.1) goto 140
        same_layer = .true.
        j = l_seq(1)
        i = 2
  150   if(l_seq(i).eq.j) then
          i = i + 1
          if(i.le.l_cnt) goto 150
        else
          same_layer = .false.
        endif
*
* Check if any of the layer transitions have non-zero probability
* initialize flags so that we do not produce duplicate error messages
        do 160 i = 1, n_layers
          do 170 j = 1, n_layers
            did_it(j,i) = .false.
  170     continue
  160   continue
*
* now check for legal transitions
        write(op,400)
        do 180 n = 1, l_cnt-1
          i = l_seq(n)
          j = l_seq(n+1)
          if(.not.there(j,i)) then
            ok = .false.
            if(.not.did_it(j,i)) then
              did_it(j,i) = .true.
              write(op,401) j, i
            endif
          endif
  180   continue
  140   continue
      endif
      if(.not.ok) goto 999
*
* Pre-compute a pseudo-Lorentzian form factor for lateral
* planar widths.
*
      if(finite_width) then
* ideally, FFACT_SIZE is odd
        m = FFACT_SIZE/2
* incr contains the correct increment size such that the first and
* last array elements contain zero.
        incr = (THREE*N_SIGMAS*N_SIGMAS + ONE) / (TWO*N_SIGMAS*m)
        ffact_scale = incr
        z = ONE + N_SIGMAS*N_SIGMAS
        tmp = ONE / (z*z)
* put the peak in the middle (if FFACT_SIZE is odd)
        formfactor(m+1) = ONE
        do 190, n = 1, m
          z = n*incr
          if(z.le.dble(N_SIGMAS)) then
* Lorentzian part
            x = ONE / (ONE + z*z)
          else
* Linear part
             x = tmp*(THREE*N_SIGMAS*N_SIGMAS + ONE - TWO*N_SIGMAS*z)
          endif
          if(m+n+1.le.FFACT_SIZE) formfactor(m+n+1) = x
          if(m-n+1.gt.0)          formfactor(m-n+1) = x
  190   continue
* compute half width in reciprocal Angstroms for use in computing the
* subtle impact of a-b shape broadening on 00l reflections
        tmp = Wa*sin(PI-cell_gamma)
        ffwdth = sqrt(ONE/(tmp*tmp) + ONE/(Wb*Wb))
      endif
*
* get diffraction symmetry
*
* establish some bounds
      no_trials = 25
      max_var = ZERO
* save lambda, HKL_LIM (called by THRESH), may change it.
      old_lambda = lambda
*
      call THRESH(ok)
      if(.not.ok) goto 999
*
      if(SymGrpNo.eq.11) then
        write(op,202)
     |      'Axial integration only selected.'
      else
        write(op,202)
     |      'Evaluating point group symmetry of diffraction data. . .'
      endif
*      if(DoSymDump) then
*        call GETFNM(rootnam, sym_fnam, 'sym', ok)
*        if(.not.ok) then
*          write(op,202) 'OPTIMZ: ERROR in creating symmetry dumpfile.'
*          goto 999
*        endif
*        if(sy.ne.op)
*     |      open(unit = sy, file = sym_fnam, status = 'new')
*        write(op,204) 'Writing symmetry data to file ''',
*     |                 sym_fnam(1:LENGTH(sym_fnam)),'''. . .'
*        write(sy,303) '''', rootnam(1:LENGTH(rootnam)),''''
*        write(sy,203) no_trials
*        write(sy,206) tiny_inty
*      endif
*
      check_sym = .false.
      if(SymGrpNo.eq.UNKNOWN) then
        SymGrpNo = GET_SYM(ok)
        if(.not.ok) goto 999
        write(op,200) 'Diffraction point symmetry is ',pnt_grp
        if(SymGrpNo.ne.1) then
          write(op,201) '  to within a tolerance of one part in ',
     |                  nint(ONE / tolerance)
        endif
      else
        check_sym = .true.
        call CHK_SYM(ok)
        if(.not.ok) goto 999
      endif
*
* restore lambda.
      lambda = old_lambda
*
      has_l_mirror = SymGrpNo.ne.1 .and. SymGrpNo.ne.3 .and.
     |               SymGrpNo.ne.5 .and. SymGrpNo.ne.6 .and.
     |               SymGrpNo.ne.11
*
      if(DoSymDump) then
        write(sy,202) ' '
        write(sy,204)
     | 'The diffraction data fits the point group symmetry ''',
     | pnt_grp(1:LENGTH(pnt_grp)),''''
        if(SymGrpNo.ne.1 .and. SymGrpNo.ne.11) then
          if(max_var.gt.eps6 .and. max_var.le.eps1) then
            write(sy,201) '  with a tolerance of one part in ',
     |                  nint(ONE / max_var)
          else if(max_var.gt.eps1) then
            write(sy,205) '  with a tolerance of one part in ',
     |                  ONE / max_var
          else
            write(sy,202)
     |      '  with a tolerance better than one part in a million.'
          endif
        else
          write(sy,202)
     |'By definition, all diffraction data has a center of symmetry'
          write(sy,202) 'thus, there is no need to test for inversion.'
        endif
* close file, unless the output was to the default device
        if(sy.ne.op) close (unit = sy)
      endif
* establish integration limits and weighting factors
      call GET_BDS()
* compute angles of scanning vectors relative to 1 0 0
      theta1 = HKANGL(k_start, h_start)
      theta2 = HKANGL(k_end, h_end)
* resolve ambiguity in the case of -1 point symmetry
      if(SymGrpNo.eq.1 .or. SymGrpNo.eq.11) then
        theta1 = -PI
        theta2 =  PI
      endif
  999 return
  200 format(1x, 2a)
  201 format(1x, a, i6)
  202 format(1x, a)
  203 format(1x, 'Number of trials per symmetry element = ', i4)
  204 format(1x, 3a)
  205 format(1x, a, f3.1)
  206 format(1x, "Threshold intensity = ", g22.6)
  302 format(1x, 'SYMMETRY EVALUATIONS FOR DATA IN FILE ''', a, '''')
  303 format(1x, 'SYMMETRY EVALUATIONS FOR DATA IN FILE ', 3a)
  400 format(1x,'Checking for conflicts in layer stackings . . .')
  401 format(1x,'ERROR: Layer ',i2,' cannot stack after layer ',i2)
  500 format(1x,'ERROR: Non-physical interlayer uncertainty factor ',a)
  501 format(1x,'       for stacking from layer ', i2, ' to ', i2)
  502 format(1x,'       It is too large relative to ', a)
      end
*
* ______________________________________________________________________
* Title: OVERLP
* Authors: MMJT
* Date: 24 Feb 1990
* Description: This function compares the coordinates of atoms in
* adjacent layers, and searches for any overlap. If it finds any
* overlap, it checks the summed occupancy to check that it does
* not exceed 1. If OVERLP finds atoms too close it provides a warning
* message only.
*
*      ARGUMENTS:
*            No arguments are passed.
*
*      COMMON VARIABLES:
*            uses:  n_layers, there, l_n_atoms, l_r, a_pos
*
*        modifies:  No COMMON variables are modified
* ______________________________________________________________________
*
      subroutine OVERLP()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical invert
      character*33 txt
      integer*4 i, j, m, n, nn, j2, err_no, max_err, fact, at_num
      integer*4 PRUNE
      parameter(max_err = 100)
      real*8 lay(3,2*MAX_A)
      real*8 x1, y1, z1, x2, y2, z2, sum_occ, tol, tmp, BOUNDS
      parameter(tol = eps1)
*
* external functions
      external BOUNDS, PRUNE
*
      write(op,400) 'Checking for conflicts in atom positions . . .'
*
      err_no = 0
      do 10 i = 1, n_layers
        fact = 1
        invert = l_symmetry(i).eq.CENTRO
        if(invert) fact = 2
        at_num = l_n_atoms(i)
        do 20 j2 = 1, at_num
          lay(1,j2) = BOUNDS(a_pos(1,j2,i))
          lay(2,j2) = BOUNDS(a_pos(2,j2,i))
* Remember, only the a and b directions are truly periodic.
* The scaling along c is arbitrary. Furthermore, c is perpendicular
* to a and b, and is not necessarily parallel to Rii, which (if it
* exists) would define the third cell-repeat direction. In other words,
* for the general case, we cannot use the BOUNDS function along c.
          lay(3,j2) = a_pos(3,j2,i)
   20   continue
        if(invert) then
          do 30 j2 = 1, at_num
            lay(1,at_num+j2) = BOUNDS(-a_pos(1,j2,i))
            lay(2,at_num+j2) = BOUNDS(-a_pos(2,j2,i))
            lay(3,at_num+j2) = -a_pos(3,j2,i)
   30     continue
        endif
        do 40 m = 1, at_num
          x1 = lay(1,m)
          y1 = lay(2,m)
          z1 = lay(3,m)
          do 50 n = m + 1, fact * at_num
            if(n.gt.at_num) then
              nn = n - at_num
            else
              nn = n
            endif
            x2 = lay(1,n)
            y2 = lay(2,n)
            z2 = lay(3,n)
            if(abs(x1-x2)*cell_a.gt.tol) goto 50
            if(abs(y1-y2)*cell_b.gt.tol) goto 50
            if(abs(z1-z2)*cell_c.le.tol) then
              sum_occ = a_occup(nn,i) + a_occup(m,i)
              if((sum_occ - ONE).gt.eps4) then
                if(n.le.at_num) then
                  txt = 'are too close in layer'
                else
                  txt = '(inverted) are too close in layer'
                endif
                write(op,410) 'Atom ', a_name(nn,i), a_number(nn,i),
     |                   ' and atom ', a_name(m,i), a_number(m,i)
                write(op,412) txt(1:PRUNE(txt)), i
                write(op,420) 'Their combined occupancy is ', sum_occ
                err_no = err_no + 1
                if(err_no.gt.max_err) goto 999
              endif
            endif
   50     continue
   40   continue
   10 continue
*
* now let's look at i-j layer transitions and generate a simple warning
* message if it seems that the layers are intertwined.
      do 60 i = 1, n_layers
        do 70 j = 1, n_layers
          if(there(j,i) .and. i.ne.j) then
            tmp = l_r(3,j,i) +
     |              low_atom(l_actual(j)) - high_atom(l_actual(i))
            if(tmp*cell_c.le.-tol) then
              write(op,430) 'Atoms from layer', j,
     |                    ' extend into layer', i
            endif
          endif
   70   continue
   60 continue
*
      return
  999 write(op,401) 'WARNING: Number of errors exceeds ', max_err
      return
  400 format(1x, a)
  401 format(1x, a, i5)
  410 format(1x, 'WARNING: ', 2(2a, ' (number ', i3, ')' ) )
  412 format(10x, a, 1x, i2)
  420 format(1x, a, g12.5)
  430 format(1x, 'WARNING: ', 2(a, i3))
      end
*
*
* ______________________________________________________________________
** Title: PARENT
** Author: MMJT
** Date: 10 Aug 1989
** Description: This routine detects whether or not a character
** string contains any parentheses, '(' and ')'. If it does, they are
** removed. PARENT is used by RDPRMS to detect whether or not the
** anisotropic temperature factors are being used.
**
**      ARGUMENTS:
**            line  -  Line of characters to parse. (input).
**                  -  Line stripped of any parentheses. (output).
**
**      PARENT returns 0 if no parentheses were found, 1 if a balanced
**      pair of parentheses were found, and -1 if there was an error.
** ______________________________________________________________________
**
*      integer*4 function PARENT(line)
*      implicit none
**
*      character*(*) line
**
*      integer*4 i, j
**
*      PARENT = 0
*      i = index(line,'(')
*      j = index(line,')')
** no parentheses
*      if(i.eq.0 .and. j.eq.0) goto 200
** right hand parenthesis occurred before the left hand parenthesis!
*      if(j.lt.i) goto 999
** only one parenthesis
*      if(i.ne.j .and. (i.eq.0 .or. j.eq.0)) goto 999
** blot out the two parentheses found
*      write(line(i:i),'(a)') ' '
*      write(line(j:j),'(a)') ' '
** check there are not any more
*      i = index(line,'(')
*      j = index(line,')')
*      if(i.eq.0 .and. j.eq.0) then
** there were no more parentheses. Data is probably ok.
*        PARENT = 1
*      else
** there were more parentheses. Data is suspect.
*        goto 999
*      endif
**
*  200 return
*  999 PARENT = -1
*      return
*      end
**
** ______________________________________________________________________
* Title: PNTINT
* Author: MMJT
* Date: 30 July 1989
* Gets the intensity at the point h, k, l. This differs from the
* subroutine POINT only in that PNTINT does not interact with the user.
*
*      ARGUMENTS:
*            h   -  reciprocal vector h-component. (input).
*            k   -  reciprocal vector k-component. (input).
*            l   -  reciprocal lattice vector l-component. (input).
*            ok  -  logical flag indicating all went well. (output).
*
*      COMMON VARIABLES:
*            uses:  a0, b0, c0, d0, lambda, recrsv, rad_type, X_RAY
*
*        modifies:  no COMMON variables are modified
*
*      PNTINT returns the intensity at h, k, l
* ______________________________________________________________________
*
      real*8 function PNTINT(h, k, l, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k
      real*8 l
      logical ok
*
      real*8 S, ANGLE, W4, INTENS, INTEN2, theta, x
      complex*16 f(MAX_L)
*
* external functions
      external INTENS, INTEN2
* external subroutines (Some compilers need them declared external)
*      external GET_F, XYPHSE, PRE_MAT
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
* ANGLE is the Bragg angle (in radians) of the h,k,l plane
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
* W4 is the X-ray polarization factor
      W4(theta) = HALF * (ONE + (cos(TWO*theta))**2)
*
      call XYPHSE(h, k)
      call PRE_MAT(h, k)
      call GET_F(f, S(h,k,l), l)
      if(recrsv) then
        x = INTENS(f, h, k, l, ok)
      else
        x = INTEN2(f, h, k, l, ok)
      endif
      if(.not.ok) then
        if(recrsv) then
          write(op,200) 'INTENS'
        else
          write(op,200) 'INTEN2'
        endif
        write(op,201) 'h,k,l = ', h,',', k,',', l
      endif
      if(rad_type.eq.X_RAY)  x = x * W4(ANGLE(h,k,l))
*
      PNTINT = x
*
      return
  200 format(1x, 'ERROR returned from ', a, ' in subroutine PNTINT')
  201 format(1x, 2(a, i3), a, g12.5)
      end
*
* ______________________________________________________________________
** Title: POINT
** Author: MWD and MMJT
** Date: 18 Aug 1988
** Description:  This subroutine prompts the user for h, k, l values
** and displays the intensity at that point.
**
**      ARGUMENTS:
**            ok  -  logical flag indicating all went well. (output).
**
**      COMMON VARIABLES:
**            uses:  a0, b0, c0, d0, lambda, cntrl, CFile, xplcit,
**                   recrsv, rad_type, n_layers, inf_thick, wavefn
**
**        modifies:  no COMMON variables are modified
** ______________________________________________________________________
**
*      subroutine POINT(ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      logical ok
**
*      logical GET_S, GET_S2
*      integer*4 h, k, i
*      character chr
*      real*8 l, INTENS, INTEN2, Shkl, x, ANGLE, SS, W4, theta, Q2
*      complex*16 f(MAX_L), s(MAX_L)
**
** external functions
*      external INTENS, INTEN2, GET_S, GET_S2
** external subroutines (Some compilers need them declared external)
**      external XYPHSE, PRE_MAT, GET_MAT, GET_F
**
** statement functions
** SS is the value of 1/d**2 at hkl
*      SS(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
** ANGLE is the Bragg angle (in radians) of the h,k,l plane
*      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(SS(h,k,l)))
** W4 is the X-ray polarization factor
*      W4(theta) = HALF * (ONE + (cos(TWO*theta))**2)
**
*    1 write(op,400) ' '
*      write(op,400) 'CALCULATING INTENSITY AT A POINT. . .'
*   10 write(op,400) 'Enter h, k, l'
*      read(cntrl,*,err=10)  h, k, l
*      if(CFile) write(op,401) h, k, l
** check angles are meaningful
*      Q2 = FOUR / (lambda**2)
*      Shkl = SS(h,k,l)
** Check that SS(h,k,l) is legal
*      if(Shkl.lt.ZERO) then
*        write(op,403) 'ERROR: In POINT(), 1/d**2 = ', Shkl
*        ok = .false.
*        goto 999
*      endif
*      if(Shkl.gt.Q2) then
*        write(op,402) h, k, l,' exceeds 180 degree scattering angle!'
*        write(op,400) 'Re-enter. . .'
*        goto 10
*      endif
** make sure we are not going to blow up at the origin
*      if(h.eq.0 .and. k.eq.0 .and. rad_type.eq.ELECTN) then
*        if(Shkl.le.eps4) then
*          write(op,400)
*     |   'Cannot integrate across the origin for electron radiation'
*          write(op,400) 'Re-enter. . .'
*          goto 10
*        endif
*      endif
** get intensity
*      call XYPHSE(h, k)
*      call PRE_MAT(h, k)
*      call GET_F(f, Shkl, l)
*      if(recrsv) then
*        x = INTENS(f, h, k, l, ok)
*        if(.not.ok) goto 999
** set up mat again to re-call GET_S (or GET_S2)
*        call PRE_MAT(h, k)
*        call GET_MAT(h, k, l)
*        if(inf_thick) then
** initialize s to -f, since mat is -(ident - T)
*          do 20 i = 1, n_layers
*            s(i) = - f(i)
*   20     continue
*          ok = GET_S(f, s, h, k, l)
*        else
** s is initialized within GET_S2, which then calls GET_S
*          ok = GET_S2(f, s, h, k, l)
*        endif
*      else
*        x = INTEN2(f, h, k, l, ok)
*      endif
*      if(.not.ok) goto 999
** for diagnostics write out the s values as well
*      if(rad_type.eq.X_RAY)  x = x * W4(ANGLE(h,k,l))
*      write(op,404)
*     |     '2theta = ', RAD2DEG * TWO * ANGLE(h,k,l), ' degrees'
*      if(Shkl .gt. ZERO) write(op,403) 'd = ', ONE / sqrt(Shkl)
*      if(Shkl .ge. ZERO) write(op,403) '1/d = ', sqrt(Shkl)
*      write(op,400) ' '
*      write(op,400) 'Layer scattering factors'
*      do 30 i = 1, n_layers
*        write(op,405) i, f(i)
*   30 continue
*      write(op,400) ' '
*      if(recrsv) then
*        write(op,400) 'Average scattered wavefunctions'
*        do 40 i = 1, n_layers
*          write(op,407) i, s(i)
*   40   continue
*      else
*        write(op,408) 'Crystal wavefunction', wavefn
*      endif
*      write(op,400) ' '
*      write(op,406) 'Intensity at ', h, k, l, ' = ', x
*      write(op,400) 'Another point? (y/n)'
*      read(cntrl,'(a)') chr
*      if(chr(1:1).eq.'y' .or. chr(1:1).eq.'Y') goto 1
**
*  999 return
*  400 format(1x, a)
*  401 format(1x, 2i3, g12.5)
*  402 format(1x, 2i3, g12.5, a)
*  403 format(1x, a, g12.5)
*  404 format(1x, a, g12.5, a)
*  405 format(1x, 'f(', i2, ') = (', g14.7, ',', g14.7, ')')
*  406 format(1x, a, 2i3, g12.5, a, g14.7)
*  407 format(1x, 'psi(', i2,') = (', g14.7, ',', g14.7, ')')
*  408 format(1x, a, ' psi = (', g14.7, ',', g14.7, ')')
*      end
**
** ______________________________________________________________________
* Title: POLINT
* Author: MMJT, adapted from Numerical Recipes Software
* Date: Copyright (C) 1985
* Description: Given arrays xa and ya, each of length n, and given
* a value x, this routine returns an interpolated estimate for y(x),
* and an error estimate dy. If P(x) is the polynomial of degree n-1,
* such that P(xa_i) = ya_i, i=1,....N, then the returned value
* y = P(x). This routine is called by APPR_F.
*
*      ARGUMENTS:
*            xa  -  Array of length n containing the x values
*                   corresponding to the n ya values. (input).
*            ya  -  Array of length n containing input y values. (input).
*            n   -  Length of the arrays entered. (input).
*            x   -  x position at which an interpolated value is
*                   required. (input).
*            y   -  interpolated y value. (output).
*            dy  -  estimate of the error in y. (output).
*            ok  -  logical flag indicating all went well. (output).
* ______________________________________________________________________
*
      subroutine POLINT(xa,ya,n,x,y,dy,ok)
      include 'DIFFaX.par'
*     save
*
      integer*4 n
      real*8 x, xa(n)
      complex*16 ya(n), dy, y
      logical ok
*
      integer*4 NMAX, i, m, ns
      parameter (NMAX = 10)
      real*8 dif, dift, ho, hp
      complex*16 c(NMAX), d(NMAX), w, den
*
      ns = 1
      dif = abs(x - xa(1))
      do 11 i = 1, n
        dift = abs(x - xa(i))
        if(dift.lt.dif) then
          ns = i
          dif = dift
        endif
        c(i) = ya(i)
        d(i) = ya(i)
   11 continue
      y = ya(ns)
      ns = ns - 1
      do 13 m = 1, n - 1
        do 12 i = 1, n - m
          ho = xa(i) - x
          hp = xa(i + m) - x
          w = c(i+1) - d(i)
          den = dcmplx(ho - hp, ZERO)
          if(den.eq.C_ZERO) goto 999
          den = w / den
          d(i) = hp * den
          c(i) = ho * den
   12   continue
        if(2*ns.lt.n-m) then
          dy = c(ns + 1)
        else
          dy = d(ns)
          ns = ns - 1
        endif
        y = y + dy
   13 continue
*
      return
  999 ok = .false.
      write(op,100) 'ERROR: Zero denominator in POLINT.'
      return
  100 format(1x, a)
      end
*
* ______________________________________________________________________
* Title: PRE_MAT
* Author: MMJT
* Date: 21 Mar 1990; 21 July 1997
* Description:  This subroutine premultiplies some of the factors
* needed to calculate the matrix mat ( = Identity - alphaij * Rij)
* at (h,k,0) which remain constant when integrating along a streak. The
* pre-multiplied factors are stored in mat1, and fatsWalla_hk.
*
*      ARGUMENTS:
*            h   -  reciprocal lattice vector h-component. (input).
*            k   -  reciprocal lattice vector k-component. (input).
*
*      COMMON VARIABLES:
*            uses:   n_layers, l_r, there, detune, a0, b0, ab0, r_B11,
*                    r_B22, r_B12, same_Bs, Bs_zero, PI2, l_alpha
*
*        modifies:   l_phi, mat1, fatsWalla_hk
* ______________________________________________________________________
*
      subroutine PRE_MAT(h, k)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k
*
      real*8 dot
      integer*4 i, j
*
* Set up matrix that represents the sequences
* For the matrix inversion routines, 'mat' and 'mat1' have to be
* in 'i,j' format rather than the quicker 'j,i' format
      do 10 i = 1, n_layers
        do 20 j = 1, n_layers
          if(there(j,i)) then
            dot = PI2*(h*l_r(1,j,i) + k*l_r(2,j,i))
            l_phi(j,i) = dcmplx(cos(dot),sin(dot))
            if(same_Bs.or.Bs_zero(j,i)) then
              mat1(i,j) =  detune(j,i) * l_alpha(j,i) * l_phi(j,i)
            else
* h-k components only. l-components are handled later by GET_MAT.
              mat1(i,j) = detune(j,i) * l_alpha(j,i) * l_phi(j,i)
     |           * exp(-QUARTER*(r_B11(j,i)*a0*h*h + r_B22(j,i)*b0*k*k)
     |           + HALF*r_B12(j,i)*ab0*h*k )
            endif
          else
            mat1(i,j) = C_ZERO
          endif
   20   continue
   10 continue
*
* Are all the uncertainty factors identical?
* Here, we compute only the h-k components.
* The l-components are handled later by GET_MAT.
      if(same_Bs) then
        if(all_Bs_zero) then
* This initialization is not actually necessary, since if we are here,
* fatsWalla_hk will not be needed by GET_MAT. However, let's be safe.
          fatsWalla_hk = ONE
        else
          fatsWalla_hk =
     |      exp(-(QUARTER*(a_B11*a0*h*h + a_B22*b0*k*k) +
     |          HALF*a_B12*ab0*h*k) )
        endif
      endif
*
      return
      end
*
* ______________________________________________________________________
* Title: PRUNE
* Author: MMJT
* Date: 4 Oct 1989
* Description: This function determines the position of the last
* non-blank character in a character string.
*
*      ARGUMENTS:
*            line     -  Line of characters to examine. (input).
*
* PRUNE returns the position of the last non-blank character, or zero
* if there was an error.
* ______________________________________________________________________
*
      integer*4 function PRUNE(line)
      include 'DIFFaX.par'
*     save
*
      character*(*) line
*
      integer*4 lin_len, i
*
      PRUNE = 0
*
      lin_len = len(line)
      i = lin_len
   10 if(i.gt.0) then
        if(line(i:i).eq.' ') then
          if(i.gt.1) then
            i = i - 1
            goto 10
          endif
        endif
      endif
*
      if(i.gt.0) PRUNE = i
*
      return
      end
*
* ______________________________________________________________________
* Title: PV
* Author: MMJT
* Date: 21st Dec 1990; 7 Mar 1995; 9 June 1998
* Description: This subroutine performs the pseudo-Voigt
* instrumental broadening. Expects the standard u, v and w and gamma
* parameters to be available in 'common'. PV does not conserve
* intensity well when the broadening width is comparable to d_theta.
* Data at the extreme ends of the spectrum are corrupted slightly.
*
*      ARGUMENTS:
*            th2_low  -  lowest 2theta angle to consider. (input).
*
*      COMMON VARIABLES:
*            uses:  pv_u, pv_v, pv_w, pv_gamma, d_theta, spec
*                   NONE, PI, RAD2DEG, blurring, th2_max
*
*        modifies:  brd_spc
* ______________________________________________________________________
*
      subroutine PV(th2_low)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 th2_low
*
      integer*4 i, j, n_low, n_high, indx
      real*8 th_rng, tn_th, c00, hk_inv, th0
      real*8 k1, k2, k3, k4, k5, pVoigt, const, tmp, speci
*
* first check the numbers
      if(pv_u.eq.ZERO .and. pv_v.eq.ZERO .and. pv_w.eq.ZERO) goto 990
      if(pv_gamma.lt.ZERO .or. pv_gamma.gt.ONE) goto 999
*
      if(th2_low.lt.ZERO .or. th2_low.ge.th2_max) then
        write(op,103) 'PV: Cut-off angle ', th2_low,
     |        ' is out of bounds. Angle reset to zero.'
        th2_low = ZERO
      endif
*
* th2_low is the angle relative to th2_min
* 2*d_theta is the angular step size
      n_low  = int(HALF*th2_low/d_theta) + 1
      n_high = int(HALF*(th2_max-th2_min)/d_theta) + 1
*
      c00 = FOUR * log(TWO)
      const = TWO * RAD2DEG * d_theta
      do 10 i = 1, n_high
        brd_spc(i) = ZERO
   10 continue
      k1 = pv_gamma*TWO/PI
      k2 = (ONE - pv_gamma)*sqrt(c00/PI)
      k3 = -c00
      th0 = HALF*th2_min
      do 20 i = n_low, n_high
* get tan((2theta)/2)
        tn_th = tan(i*d_theta + th0)
        tmp = (pv_u * tn_th + pv_v) * tn_th  +  pv_w
        if(tmp.le.ZERO) goto 995
        hk_inv = ONE / sqrt(tmp)
        tmp = (const * hk_inv)**2
        k4 = k1 * hk_inv * const
        k5 = k2 * hk_inv * const
        speci = spec(i)
        do 30 j = n_low - i, n_high - i
          th_rng = tmp * j*j
          pVoigt = (k4/(ONE+FOUR*th_rng) + k5*exp(k3*th_rng)) * speci
          indx = i + j
          brd_spc(indx) = brd_spc(indx) + pVoigt
   30   continue
   20 continue
      return
  990 write(op,100) 'pseudo-Voigt parameters are zero in PV()'
      write(op,101)
      blurring = NONE
      return
  995 write(op,102)
     | 'ERROR: pseudo-Voigt spread function is complex at theta = ',
     |  i*d_theta
      write(op,100) '   u, v, w parameters are illegal.'
      write(op,101)
      blurring = NONE
      return
  999 write(op,100) 'Illegal pseudo-Voigt gamma value in PV()'
      write(op,101)
      blurring = NONE
      return
  100 format(1x, a)
  101 format(1x, '   pseudo-Voigt instrumental broadening not added')
  102 format(1x, a, g12.5)
  103 format(1x, a, g12.5, a)
      end
*
* ______________________________________________________________________
** Title: RDBLUR
** Authors: MMJT
** Date: 4 Oct 1989
** Description: This function reads the type of instrumental
** broadening to be simulated, and the associated parameters.
** Legal types are 'NONE ', 'PSEUDO-VOIGT ', 'GAUSSIAN '
** and 'LORENTZIAN '. PSEUDO-VOIGT requires 4 parameters,
** GAUSSIAN requires a standard deviation in 2theta degrees, and
** LORENTZIAN requires a halfwidth in 2theta degrees.
** If the keyword 'TRIM' is at the end of the data line, intensity
** near the origin will be ignored when the instrumental broadening
** is added.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**        modifies:  blurring, pv_u, pv_v, pv_w, pv_gamma,
**                   FWHM, GAUSS, LORENZ, NONE, PS_VGT, trim_origin
**
**      RDBLUR returns logical .true. if a legal function for
**      simulating instrumental broadening was entered.
**
**      Legal types are:
**             'NONE'
**             'PSEUDO-VOIGT' (with u, v, w and gamma parameters)
**             'GAUSSIAN' (with a standard deviation in degrees)
**             'LORENTZIAN' (with a half width in degrees)
** ______________________________________________________________________
**
*      logical function RDBLUR(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line, dummy
*      character*80 list(4)
*      integer*4 CHOICE, LENGTH, CNTARG, i, iflag, arg_num, lin_len
**
** external functions
*      external CHOICE, LENGTH, CNTARG
** external subroutines (Some compilers need them declared external)
**      external GETLNE, WRTLNE, TOUPPR
**
*      RDBLUR = .false.
**
**
*      messge = 'Instrumental broadening parameters were not specified.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
*      lin_len = len(line)
**
** First, check whether or not user wants to trim the low angle
** intensities for the broadened spectrum. If present, remove the
** keyword 'TRIM', and set flag.
*      i = index(line, ' TRIM')
*      if(i.eq.0) then
*        trim_origin = .false.
*      else
*        trim_origin = .true.
*        write(line(i+1:i+4), '(a)') '    '
*      endif
**
*      list(1) = 'NONE '
*      list(2) = 'GAUSSIAN '
*      list(3) = 'LORENTZIAN '
*      list(4) = 'PSEUDO-VOIGT '
*      iflag   =  CHOICE(line, list, 4)
**
*      if(iflag.lt.1 .or. iflag.gt.4) then
*        write(op,100)
*     |  'Illegal instrumental broadening type. Legal types are:'
*        do 10 i = 1, 4
*          write(op,101) '               ', list(i)
*   10   continue
*        messge = 'Instrumental broadening incorrectly specified.$'
*        goto 999
*      endif
**
*      if(iflag.eq.1) blurring = NONE
**
*      if(iflag.eq.2) then
**
** check the number of arguments
*        dummy = line(LENGTH('GAUSSIAN')+2:lin_len)
*        arg_num = CNTARG(dummy)
*        if(arg_num.ne.1 .and. arg_num.ne.3) then
*          messge =
*     |       'Illegal number of GAUSSIAN arguments. FWHM, or u,v,w.$'
*          goto 999
*        endif
**
*        if(arg_num.eq.1) then
*          blurring = GAUSS
**
*          messge = 'Problems reading Gaussian FWHM.$'
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) dummy
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) FWHM
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**          read(dummy, *, err=999) FWHM
**
*          if(FWHM.lt.ZERO) then
*            messge = 'Gaussian FWHM is -ve.$'
*            goto 999
*          endif
**
*          if(FWHM.lt.eps7) then
*            write(op,102) 'Gaussian FWHM is zero.'
*            write(op,102)
*     |        ' No instrumental broadening will be simulated.'
*            blurring = NONE
*          endif
**
*        else if(arg_num.eq.3) then
*          blurring = PV_GSS
**
*          messge = 'Problems reading Gaussian u, v, and w factors.$'
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) dummy
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) pv_u, pv_v, pv_w
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**          read(dummy, *, err=999) pv_u, pv_v, pv_w
**
*          pv_gamma = ZERO
**
*          if(pv_w.lt.ZERO) then
*            write(op,102)
*     |        'WARNING: Gaussian w-factor is less than zero.'
*            write(op,102)
*     |        '  This may cause the broadening function to fail.'
*          endif
**
*          if(abs(pv_u).lt.eps7 .and. abs(pv_v).lt.eps7 .and.
*     |       abs(pv_w).lt.eps7) then
*            write(op,102) 'Gaussian u, v, w factors are zero.'
*            write(op,102) 'No instrumental broadening will be applied.'
*            blurring = NONE
*          endif
**
** are the parameters equivalent to a Gaussian with
** constant standard deviation?
*          if(abs(pv_u).le.eps7 .and. abs(pv_v).le.eps7 .and.
*     |                                    pv_w.gt.ZERO) then
*            FWHM = sqrt(pv_w)
*            write(op,102)
*     |       'Gaussian factors imply a constant FWHM ', FWHM
*            blurring = GAUSS
*          endif
**
*        endif
**
*      endif
**
*      if(iflag.eq.3) then
**
** check the number of arguments
*        dummy = line(LENGTH('LORENTZIAN')+2:lin_len)
*        arg_num = CNTARG(dummy)
*        if(arg_num.ne.1 .and. arg_num.ne.3) then
*          messge =
*     |     'Illegal number of LORENTZIAN arguments. FWHM, or u,v,w.$'
*          goto 999
*        endif
**
*        if(arg_num.eq.1) then
*          blurring = LORENZ
**
*          messge = 'Problems reading Lorentzian width.$'
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) dummy
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) FWHM
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**          read(dummy, *, err=999) FWHM
**
*          if(FWHM.lt.ZERO) then
*            messge =
*     |        'Negative value for Lorentzian FWHM. Must be positive.$'
*            goto 999
*          endif
**
*          if(FWHM.lt.eps7) then
*            write(op,102) 'Lorentzian FWHM is zero.'
*            write(op,102)
*     |      ' No instrumental broadening will be simulated.'
*            blurring = NONE
*          endif
**
*        else if(arg_num.eq.3) then
*          blurring = PV_LRN
**
*          messge = 'Problems reading Lorentzian u, v, and w factors.$'
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) dummy
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) pv_u, pv_v, pv_w
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**          read(dummy, *, err=999) pv_u, pv_v, pv_w
**
*          pv_gamma = ONE
**
*          if(pv_w.lt.ZERO) then
*            write(op,102)
*     |        'WARNING: Lorentzian w-factor is less than zero.'
*            write(op,102)
*     |        '  This may cause the broadening function to fail.'
*          endif
**
*          if(pv_u.lt.eps7 .and. pv_v.lt.eps7 .and.
*     |                         abs(pv_w).lt.eps7) then
*            write(op,102) 'Lorentzian u, v, w factors are zero.'
*            write(op,102) 'No instrumental broadening will be applied.'
*            blurring = NONE
*          endif
**
** are the parameters equivalent to a Gaussian with
** constant standard deviation?
*          if(pv_u.le.eps7 .and. pv_v.le.eps7 .and. pv_w.gt.ZERO) then
*            FWHM = sqrt(pv_w)
*            write(op,103)
*     |       'Lorentzian factors imply a constant FWHM ', FWHM
*            blurring = LORENZ
*          endif
**
*        endif
**
*      endif
**
*  500 continue
**
*      if(iflag.eq.4) then
*        blurring = PS_VGT
**
** check the number of arguments
*        dummy = line(LENGTH('PSEUDO-VOIGT')+2:len(line))
*        arg_num = CNTARG(dummy)
*        if(arg_num.ne.4) then
*          messge = 'Exactly 4 pseudo-Voigt arguments are required.$'
*          goto 999
*        endif
**
*        messge =
*     |      'Problems reading pseudo-Voigt u, v, w, gamma factors.$'
**
** for strict FORTRAN77 compliance
*        write(scrtch, *, err=990) dummy
*        rewind(scrtch, err=990)
*        read(scrtch, *, err=999) pv_u, pv_v, pv_w, pv_gamma
*        rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**        read(dummy, *, err=999) pv_u, pv_v, pv_w, pv_gamma
**
*        if(pv_w.lt.ZERO) then
*          write(op,102)
*     |        'WARNING: pseudo-Voigt w factor is less than zero.'
*          write(op,102)
*     |        '  This may cause the pseudo-Voigt function to fail.'
*        endif
**
*        if(abs(pv_u).lt.eps7 .and. abs(pv_v).lt.eps7 .and.
*     |                         abs(pv_w).lt.eps7) then
*          write(op,102) 'Pseudo-Voigt u, v, w factors are zero.'
*          write(op,102) 'No instrumental broadening will be simulated.'
*          blurring = NONE
*        endif
**
*        if(pv_gamma.lt.ZERO) then
*          messge =
*     |   'Pseudo-Voigt gamma factor is -ve (must be between 0 and 1).$'
*          goto 999
*        endif
**
*        if(pv_gamma.gt.ONE) then
*          messge = 'Pseudo-Voigt gamma factor is greater than 1.$'
*          goto 999
*        endif
**
** are the pseudo-Voigt parameters equivalent to a Gaussian with
** constant standard deviation?
*        if(abs(pv_u).le.eps7 .and. abs(pv_v).le.eps7 .and.
*     |     pv_w.gt.ZERO .and. pv_gamma.le.eps4) then
*          FWHM = sqrt(pv_w)
*          write(op,102)
*     |        'Pseudo-Voigt factors are equivalent to a Gaussian'
*          write(op,103) ' with FWHM ', FWHM
*          blurring = GAUSS
*        endif
**
** are the pseudo-Voigt parameters equivalent to a Lorentzian with
** constant half-width?
*        if(abs(pv_u).le.eps7 .and. abs(pv_v).le.eps7 .and.
*     |     pv_w.gt.ZERO .and. abs(pv_gamma - ONE).le.eps4) then
*          FWHM = sqrt(pv_w)
*          write(op,102)
*     |        'Pseudo-Voigt factors are equivalent to a Lorentzian'
*          write(op,103) ' with FWHM ', FWHM
*          blurring = LORENZ
*        endif
**
*      endif
**
*      RDBLUR = .true.
*      return
*  990 messge = 'Problems using scratch file in RDBLUR.$'
*  999 call WRTLNE(line)
*      return
*  100 format(1x, 'ERROR: ', a)
*  101 format(1x, 2a)
*  102 format(1x, a)
*  103 format(1x, a, g12.5)
*      end
**
** ______________________________________________________________________
** Title: RDCELL
** Authors: MMJT
** Date: 4 Oct 1989
** Description: This function reads the layer mesh dimensions from
** the second data line of the data file. It checks that the header
** 'STRUCTURAL' is present.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**        modifies:  cell_a, cell_b, cell_c, cell_gamma
**
**      RDCELL returns logical .true. if acceptable values
**      for the layer mesh dimensions were read.
** ______________________________________________________________________
**
*      logical function RDCELL(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line
**
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR, WRTLNE
**
*      RDCELL = .false.
**
** first check that 'STRUCTURAL' header is present
*      messge = '''STRUCTURAL'' header could not be read.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
*      messge = '''STRUCTURAL'' header not found.$'
*      if(index(line,'STRUCTURAL ').eq.0) goto 999
**
*      messge = 'a, b, c, and gamma are improperly specified.$'
*      call GETLNE(unit_no, line, *999)
**
** for strict FORTRAN77 compliance
*      write(scrtch, *, err=990) line
*      rewind(scrtch, err=990)
*      read(scrtch, *, err=999) cell_a, cell_b, cell_c, cell_gamma
*      rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**      read(line, *, err=999) cell_a, cell_b, cell_c, cell_gamma
**
*      if(cell_a.le.ZERO) then
*        messge = 'Negative (or zero) value for cell_a.$'
*        goto 999
*      endif
**
*      if(cell_b.le.ZERO) then
*        messge = 'Negative (or zero) value for cell_b.$'
*        goto 999
*      endif
**
*      if(cell_c.le.ZERO) then
*        messge = 'Negative (or zero) value for cell_c.$'
*        goto 999
*      endif
**
*      if(cell_gamma.le.ZERO) then
*        messge = 'Negative (or zero) value for cell_gamma.$'
*        goto 999
*      endif
**
*      if(cell_gamma.ge.ONE_EIGHTY) then
*        messge = 'Cell_gamma is greater than 180 degrees.$'
*        goto 999
*      endif
**
*      cell_gamma = cell_gamma * DEG2RAD
*      RDCELL = .true.
*      return
*  990 messge = 'Problems using scratch file in RDCELL.$'
*  999 call WRTLNE(line)
*      return
*      end
**
** ______________________________________________________________________
** Title: RDFILE
** Author: MMJT
** Date: 3 Oct 1989; 3 Mar 1995
** Description: This function reads in the data from the user-defined
** data file. The data file is in Naur-Backus form, and the user may
** insert '{' and '}' delimiters to comment the data. The data file
** is not read in its raw form. The subroutine GETLNE is used
** to read each line, stripped of its comments, and left-justified.
**
**      ARGUMENTS:
**           fname  -  The name of the input data file to read. (input).
**
**      RDFILE returns logical .true. if the file read went well.
** ______________________________________________________________________
**
*      logical function RDFILE(fname)
*      include 'DIFFaX.par'
**     save
**
*      character*(*) fname
**
*      character*200 messge
*      logical goodfile
*      logical RDRADN, RDWAVL, RDBLUR, RDCELL, RDPTGP
*      logical RDNLAY, RDWDTH, RDLAYR, RDSTAK, RDTRNS
*      integer*4 LENGTH, unit_no
**
** external functions
*      external RDRADN, RDWAVL, RDBLUR, RDCELL, RDPTGP
*      external RDNLAY, RDWDTH, RDLAYR, RDSTAK, RDTRNS
*      external LENGTH
**
*      RDFILE = .false.
*      unit_no = df
*      write(op,100) 'Reading ''', fname(1:LENGTH(fname)),''''
**
** now read the file, one datum at a time.
*                   goodfile = RDRADN(unit_no, messge)
*      if(goodfile) goodfile = RDWAVL(unit_no, messge)
*      if(goodfile) goodfile = RDBLUR(unit_no, messge)
*      if(goodfile) goodfile = RDCELL(unit_no, messge)
*      if(goodfile) goodfile = RDPTGP(unit_no, messge)
*      if(goodfile) goodfile = RDNLAY(unit_no, messge)
*      if(goodfile) goodfile = RDWDTH(unit_no, messge)
*      if(goodfile) goodfile = RDLAYR(unit_no, messge)
*      if(goodfile) goodfile = RDSTAK(unit_no, messge)
*      if(goodfile) goodfile = RDTRNS(unit_no, messge)
**
*      if(goodfile) then
*        write(op,100) '''', fname(1:LENGTH(fname)),''' read in.'
*      else
*        write(op,200) messge(1:index(messge,'$') - 1)
*      endif
**
*      write(messge,101)
*     |      'Problems closing ''', fname(1:LENGTH(fname)),'''$'
*      close(df,err=900)
** if all went well, set RDFILE to .true.
*      if(goodfile) RDFILE = .true.
*      return
**
*  900 write(op,200) messge(1:index(messge,'$') - 1)
*      return
*  100 format(1x, 3a)
*  101 format(3a)
*  200 format(1x, 'ERROR: ', a)
*      end
**
** ______________________________________________________________________
** Title: RDLAYR
** Authors: MMJT
** Date: 7 June 1990
** Description: This function reads in the layer atomic coordinates.
** Each set of data must be prefaced with a layer number, and a
** statement of whether or not the layer is centrosymmetric. If the
** layer is centrosymmetric, only the asymmetric coordinates need
** be entered (with appropriate fractional occupancies for special
** positions.
** If a layer has the same structure as a previously defined layer
** the layer can be defined as 'LAYER i = j', where j is a previously
** defined layer.
** The format assumed for atom coordinates is:
**
**    'atom name', number, x, y, z, Debye-Waller factor, occupancy
**
** 'atom name' must occupy 4 characters, eg. 'Si4+', 'O 1-' etc...
** The allowed atom names can be found in the file 'data.sfc'.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**            uses:  n_layers, l_actual, CENTRO, NONE
**
**        modifies:  l_symmetry, a_name, a_number, a_pos, a_B, a_occup,
**                   l_n_atoms, n_actual
**
**      RDLAYR returns logical .true. if all the layer data
**      was read successfully.
** ______________________________________________________________________
**
*      logical function RDLAYR(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line, tmpline, arg
*      character*80 list(2)
*      integer*4 CHOICE, LENGTH, NXTARG, m, n, i, j, i2, iflag
*      real*8 RDNMBR, tmp
*      logical*4 ok
**
** external functions
*      external CHOICE, LENGTH, NXTARG, RDNMBR
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR, WRTLNE
**
*      RDLAYR = .false.
**
*      do 5 i = 1, n_layers
*        high_atom(i) = ZERO
*        low_atom(i)  = ZERO
*    5 continue
**
*      m = 0
*      messge = 'Unexpected EOF before LAYER 1.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
**
*      do 30 i = 1, n_layers
*        write(messge,100) 'LAYER header not seen for LAYER ',i,'.$'
*        if(index(line, 'LAYER').ne.1) goto 999
** read layer number
*        write(messge,100) 'LAYER ',i,' number sequence invalid.$'
*        write(tmpline, '(a)', err=999)
*     |       line(LENGTH('LAYER')+1:len(line))
** for strict FORTRAN77 compliance
*        write(scrtch, *, err=990) tmpline
*        rewind(scrtch, err=990)
*        read(scrtch, *, err=999) j
*        rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**        read(tmpline, *, err=999) j
**
** layers must be entered in sequence
*        if(j.ne.i) goto 999
**
*        j = index(line, '=')
** see if this layer is equivalenced to another
*        if(j.ne.0) then
**
*          write(messge,100)
*     |           'Invalid equivalence number for LAYER ',i,'.$'
*          write(tmpline, '(a)', err=999) line(j+1:len(line))
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) tmpline
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) i2
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**        read(tmpline, *, err=999) i2
**
** layer i2 must be defined already
*          if(i2.ge.i .or. i2.lt.1) goto 999
*          l_actual(i) = l_actual(i2)
**
*          write(messge,100) 'Unexpected EOF in LAYER ',i,'.$'
*          call GETLNE(unit_no, line, *999)
*          call TOUPPR(line)
**
*        else
*          m = m + 1
*          l_actual(i) = m
**
*          write(messge,100)
*     |           'No symmetry statement for LAYER ',i,'.$'
*          call GETLNE(unit_no, line, *999)
*          call TOUPPR(line)
**
*          list(1) = 'NONE '
*          list(2) = 'CENTROSYMMETRIC '
*          iflag = CHOICE(line, list, 2)
**
*          write(messge,100)
*     |           'Invalid symmetry statement for LAYER ',i,'.$'
*          if(iflag.lt.1 .or. iflag.gt.2) goto 999
**
*          if(iflag.eq.1) then
*            l_symmetry(m) = NONE
*          else
*            l_symmetry(m) = CENTRO
*          endif
**
*          write(messge,100) 'No atoms specified for LAYER ',i,'.$'
*          call GETLNE(unit_no, line, *999)
*          if(index(line,'LAYER').eq.1) goto 999
**
*          j = 0
*   20       j = j + 1
*            if(j.gt.MAX_A) then
*              write(op,101) 'ERROR: Too many atoms on LAYER ',i,'.'
*              write(op,102) '       Maximum allowed = ', MAX_A
*              messge = ' $'
*              goto 999
*            endif
**
*            a_name(j,m) = line(1:4)
**
*            write(messge,103)
*     |        'Problems reading data for atom ''', a_name(j,m), '''.$'
*            write(tmpline, '(a)', err=999) line(5:len(line))
** read atom identifier
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            read(arg,'(i3)') a_number(j,m)
**
** read x_relative position
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            a_pos(1,j,m) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**
** read y_relative position
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            a_pos(2,j,m) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**
** read z_relative position
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            a_pos(3,j,m) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**
** read Debye-Waller factor
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            a_B(j,m) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**
** read occupancy
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            a_occup(j,m) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**
** Check values
*            if(a_B(j,m).lt.ZERO) then
*              messge = 'Negative Debye-Waller factor.$'
*              goto 999
*            endif
**
*            if(a_occup(j,m).lt.ZERO) then
*              messge = 'Negative atom occupancy.$'
*              goto 999
*            endif
**
*            if(a_occup(j,m).gt.ONE) then
*              messge = 'Occupancy greater than 1.$'
*              goto 999
*            endif
**
** get extreme upper and lower atom positions for error-checking later on
*            tmp = a_pos(3,j,m)
*            if(tmp.gt.high_atom(m)) high_atom(m) = tmp
*            if(tmp.lt.low_atom(m))   low_atom(m) = tmp
**
*            write(messge,103)'Unexpected EOF after atom ',line(1:4),'$'
*            call GETLNE(unit_no, line, *999)
** see if all atomic specifications for this layer have been read in
** but first we must convert to uppercase for the test. If there are
** more atoms, then we must convert 'line' to its original format.
*            tmpline = line
*            call TOUPPR(line)
** check data is in correct sequence
*            messge = 'STACKING data is missing.$'
*            if(index(line,'PARAMETERS').eq.1) goto 999
*            if(index(line,'LAYER').ne.1
*     |             .and. index(line,'STACKING').ne.1) then
*              line = tmpline
*              goto 20
*            endif
*            l_n_atoms(m) = j
*        endif
*   30 continue
*      n_actual = m
**
** make sure we have genuine lowest and highest atoms in each layer.
*      do 40 i = 1, n_actual
*        if(l_symmetry(i).eq.CENTRO) then
*          if(-low_atom(i).gt.high_atom(i)) then
*            high_atom(i) = -low_atom(i)
*          else if(-low_atom(i).lt.high_atom(i)) then
*            low_atom(i) = -high_atom(i)
*          endif
*        endif
*   40 continue
**
** If we have reached here, then we have read one line too many.
*      write(messge,100) 'Problems ''backspacing'' unit ',unit_no, '.$'
*      backspace(unit = unit_no, err = 980)
**
*      RDLAYR = .true.
*  980 return
*  990 messge = 'Problems using scratch file in RDLAYR.$'
*  999 call WRTLNE(line)
*      return
*  100 format(a, i2, a)
*  101 format(1x, a, i2, a)
*  102 format(1x, a, i4)
*  103 format(3a)
*      end
**
** ______________________________________________________________________
** Title: RDNLAY
** Authors: MMJT
** Date: 4 Oct 1989
** Description: This function reads the first line of the data file
** to extract the number of layers in the structure.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**        modifies:  n_layers
**
**      RDNLAY returns logical .true. if an acceptable value
**      for the number of layers was read.
** ______________________________________________________________________
**
*      logical function RDNLAY(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line
**
** external subroutines (Some compilers need them declared external)
**      external GETLNE, WRTLNE
**
*      RDNLAY = .false.
**
*      messge = 'No lines in file.$'
*      call GETLNE(unit_no, line, *999)
**
*      messge = 'Problems reading number of layers.$'
** for strict FORTRAN77 compliance
*      write(scrtch, *, err=990) line
*      rewind(scrtch, err=990)
*      read(scrtch, *, err=999) n_layers
*      rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**      read(line, *, err=999) n_layers
**
*      if(n_layers.eq.0) then
*        messge = 'Zero number of layers entered.$'
*        goto 999
*      else if(n_layers.lt.0) then
*        messge = 'Negative number of layers entered.$'
*        goto 999
*      else if(n_layers.gt.MAX_L) then
*        write(messge,100)
*     |  'Too many layers: number should not exceed ',MAX_L,'.$'
*        goto 999
*      endif
**
*      RDNLAY = .true.
*      return
*  990 messge = 'Problems using scratch file in RDNLAY.$'
*  999 call WRTLNE(line)
*      return
*  100 format(a, i2, a)
*      end
**
** ______________________________________________________________________
** Title: RDNMBR
** Author: MMJT
** Date: 21 January 2005
** Description: Reads the character string 'numberstring' to extract a
** floating point number. The character string can be either a true
** number, such as '0.3333', or can be expressed as a ratio, '1/3'.
**
**      ARGUMENTS:
**            arg   -  Character string containing numeric data.(input).
**            ok    -  True of all went well.(output).
**
**      Returns the floating point number
** ______________________________________________________________________
**
*      real*8 function RDNMBR(numberstring,ok)
*      include 'DIFFaX.par'
*      save
**
*      character*(*) numberstring
*      logical ok
**
*      character*200 arg
*      integer*4 numerator, denominator, i, j, lin_len
*      integer*4 NXTARG
**
** external functions
*      external NXTARG
**
*      ok = .true.
*      RDNMBR = ZERO
**
*      lin_len = len(numberstring)
*      j = index(numberstring, '/')
*      if(j.eq.0) then
** It was already a number
*        read(numberstring, *) RDNMBR
*      else if(j.lt.0 .or. j.ge.lin_len) then
** An error happened
*        goto 999
*      else
** It was expressed as a ratio. Blank out the slash, '/'.
*        numberstring(j:j) = ' '
** 'numberstring' now contains two arguments, numerator and denominator.
*        i = NXTARG(numberstring,arg)
*        if(i.le.0) goto 999
*        read(arg, *) numerator
*        i = NXTARG(numberstring,arg)
*        if(i.le.0) goto 999
*        read(arg, *) denominator
*        if(denominator.le.ZERO) goto 999
*        RDNMBR = dble(numerator) / dble(denominator)
*      endif
**
*  100 return
*  999 ok = .false.
*      return
*      end
**
** ______________________________________________________________________
** Title: RDPTGP
** Author: MMJT
** Date: 4 Oct 1989
** Description: This function reads the user's estimate of the
** point group symmetry of the diffraction data.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**        modifies:  pnt_grp, SymGrpNo, tolerance
**
**      RDPTGP returns logical .true. if a legal point group
**      was entered, or if 'unknown' (with a search tolerance)
**      was specified.
**
**      Legal entries are:
**                     '-1'
**                     '2/M(1)'    (diad along streaks)
**                     '2/M(2)'    (mirror contains streaks)
**                     'MMM'
**                     '-3'
**                     '-3M'
**                     '4/M'
**                     '4/MMM'
**                     '6/M'
**                     '6/MMM'
**                     'AXIAL'     (for integration along 00l only)
**                     'UNKNOWN'   (followed by an optional tolerance)
** ______________________________________________________________________
**
*      logical function RDPTGP(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line, oldline
*      character*80 list(12)
*      integer*4 CHOICE, LENGTH, PRUNE, i, iflag
**
** external functions
*      external CHOICE, LENGTH, PRUNE
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR, WRTLNE
**
*      RDPTGP = .false.
**
*      messge = 'Illegal diffraction point group symmetry.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
**
*      list(1)  = '-1 '
*      list(2)  = '2/M(1) '
*      list(3)  = '2/M(2) '
*      list(4)  = 'MMM '
*      list(5)  = '-3 '
*      list(6)  = '-3M '
*      list(7)  = '4/M '
*      list(8)  = '4/MMM '
*      list(9)  = '6/M '
*      list(10) = '6/MMM '
*      list(11) = 'AXIAL '
*      list(12) = 'UNKNOWN '
*      iflag = CHOICE(line, list, 12)
*      if(iflag.lt.1 .or. iflag.gt.12) then
*        write(op,100) 'Unrecognized symmetry type. Legal types are:'
*        do 10 i = 1, 12
*          write(op,101) '               ',list(i)
*   10   continue
*        goto 999
*      endif
**
*      pnt_grp = list(iflag)
**
** set default tolerance on intensity variability
*      tolerance = ONE
**
*      if(iflag.ge.1 .and. iflag.le.11) then
*        SymGrpNo = iflag
*      else if(iflag.eq.12) then
** a value of UNKNOWN = -1 alerts OPTIMZ that user entered 'UNKNOWN'
*        SymGrpNo = UNKNOWN
** Strip off the symmetry point group text string.
*        oldline = line
*        write(line,'(a)') oldline(LENGTH(list(iflag))+1:len(oldline))
** If the remainder of the line is blank, no tolerance was specified.
*        i = PRUNE(line)
*        if(i.le.0) then
*          messge = 'Problems reading symmetry tolerance parameter.$'
*          goto 999
*        else if(i.gt.1) then
** There was something on the line, try to read tolerance
*          messge = 'Problems reading symmetry tolerance parameter.$'
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) line
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) tolerance
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**          read(line, *, err=999) tolerance
**
*          if(tolerance.le.ZERO) then
*            messge = 'Negative (or zero) value for tolerance.$'
*            goto 999
*          endif
**
*          if(tolerance.lt.eps2) then
*            tolerance = eps2
*            write(op,102)
*            write(op,103) tolerance
*          endif
**
*        endif
*      endif
**
** convert from a percentage
*      tolerance = tolerance * eps2
**
*      RDPTGP = .true.
*      return
*  990 messge = 'Problems using scratch file in RDPTGP.$'
*  999 call WRTLNE(line)
*      return
*  100 format(1x, 'ERROR: ', a)
*  101 format(1x, 2a)
*  102 format(1x, 'WARNING: Tolerance for symmetry tests is too small.')
*  103 format(1x, '         Tolerance is reset to ', g12.5, ' percent')
*      end
**
** ______________________________________________________________________
** Title: RDRADN
** Authors: MMJT
** Date: 15 Feb 1990
** Description: This function reads the radiation type. The choices are
** X-RAY, NEUTRON and ELECTRON. Checks that the keword 'INSTRUMENTAL'
** is present.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**            uses:  ELECTN, NEUTRN, X_RAY
**
**        modifies:  rad_type
**
**      RDRADN returns logical .true. if an acceptable
**      radiation type was read.
** ______________________________________________________________________
**
*      logical function RDRADN(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line
*      character*80 list(3)
*      integer*4 CHOICE, iflag
**
** external function
*      external CHOICE
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR, WRTLNE
**
*      RDRADN = .false.
**
** first check that 'INSTRUMENTAL' header is present
*      messge = '''INSTRUMENTAL'' header could not be read.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
*      messge = '''INSTRUMENTAL'' header not found.$'
*      iflag  =index(line,'INSTRUMENTAL ')
*      if(iflag.eq.0) goto 999
**
*      messge = 'radiation type not specified.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
*      list(1) = 'X-RAY '
*      list(2) = 'NEUTRON '
*      list(3) = 'ELECTRON '
*      iflag = CHOICE(line, list, 3)
**
*      if(iflag.eq.1) then
*        rad_type = X_RAY
*      else if(iflag.eq.2) then
*        rad_type = NEUTRN
*      else if(iflag.eq.3) then
*        rad_type = ELECTN
*      else
*        messge ='Invalid radiation type (X-RAY, NEUTRON or ELECTRON).$'
*        goto 999
*      endif
**
*      RDRADN = .true.
*      return
*  999 call WRTLNE(line)
*      return
*      end
**
** ______________________________________________________________________
** Title: RDRNGE
** Authors: MMJT
** Date: 15 Feb 1990
** Description: This function reads the angular range for the spectrum.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**
**      COMMON VARIABLES:
**            uses:  DEG2RAD
**
**        modifies:  th2_min, th2_max, d_theta
**
**      RDRNGE returns logical .true. if an acceptable angular
**      range was read.
** ______________________________________________________________________
**
*      logical function RDRNGE()
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
** external subroutines (Some compilers need them declared external)
**      external GETLNE, WRTLNE
**
*      RDRNGE = .false.
**
*  123 write(op,100) 'Enter angular range:'
*      write(op,100) '  2theta min, 2theta max, 2theta increment.'
*      read(cntrl,*,err=123) th2_min, th2_max, d_theta
*      if(CFile) write(op,105) th2_min, th2_max, d_theta
**
*      if(th2_min.lt.ZERO) then
*        write(op,110) '2theta min is negative.'
*        goto 999
*      endif
**
*      if(th2_min.ge.ONE_EIGHTY) then
*        write(op,110) '2theta min exceeds 180 degrees.'
*        goto 999
*      endif
**
*      if(th2_max.le.ZERO) then
*        write(op,110) 'Negative (or zero) value for 2theta min.'
*        goto 999
*      endif
**
*      if(th2_max.gt.ONE_EIGHTY) then
*        write(op,110) '2theta max exceeds 180 degrees.'
*        goto 999
*      endif
**
** if th2_max = 180, reduce it slightly to keep us out of trouble
*      if(th2_max.eq.ONE_EIGHTY) th2_max = th2_max - eps4
**
*      if(d_theta.le.ZERO) then
*        write(op,110) 'Negative (or zero) value for 2theta increment.'
*        goto 999
*      endif
**
*      if(th2_min.ge.th2_max) then
*        write(op,110) '2theta min is larger than 2theta max.'
*        goto 999
*      endif
**
*      if((th2_max - th2_min).lt.d_theta) then
*        write(op,110)
*     |  '2theta increment is larger than 2theta max - 2theta min.'
*        goto 999
*      endif
**
*      th2_min = th2_min * DEG2RAD
*      th2_max = th2_max * DEG2RAD
*      d_theta = HALF * DEG2RAD * d_theta
**
*      RDRNGE = .true.
*      return
**
*  999 return
*  100 format(1x, a)
*  105 format(1x, g11.4, 2(3x, g11.4))
*  110 format(1x, 'ERROR: ', a)
*      end
**
** ______________________________________________________________________
** Title: RDSTAK
** Author: MMJT
** Date: 4 Oct 1989
** Description: This function reads the STACKING input of the datafile.
** The stacking is either 'EXPLICIT', whereupon DIFFaX expects the
** user to enter an explicit sequence of layers, or 'RECURSIVE'. If
** 'EXPLICIT' was specified, the user may state that the sequence is
** 'RANDOM', whereby DIFFaX will generate an explicit sequence of
** layers weighted by the stacking probabilities specified later on
** under 'PARAMETERS'. The user can specify a crystal thickness (in
** terms of the number of layers) under the 'RECURSIVE' option.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**            uses:  n_layers
**
**        modifies:  recrsv, xplcit, rndm, l_cnt, l_seq, inf_thick
**
**      RDSTAK returns logical .true. if the stacking
**      method was properly specified.
** ______________________________________________________________________
**
*      logical function RDSTAK(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line, tmpline
*      character*80 list(2)
*      integer*4 CHOICE, LAYCNT, LENGTH, i, j, iflag
**
** external functions
*      external CHOICE, LAYCNT, LENGTH
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR, WRTLNE
**
*      RDSTAK = .false.
**
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
*      if(index(line,'STACKING ').ne.1) then
*        write(messge,103) 'The ', n_layers,
*     |   ' layers are read in. ''STACKING'' section not found.$'
*        goto 999
*      endif
**
*      messge = 'Illegal STACKING description.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
**
*      list(1) = 'EXPLICIT '
*      list(2) = 'RECURSIVE '
*      iflag = CHOICE(line, list, 2)
**
*      if(iflag.lt.1 .or. iflag.gt.2) then
*        write(op,100) 'Illegal STACKING keyword. Legal types are:'
*        do 10 i = 1, 2
*          write(op,101) '               ',list(i)
*   10   continue
*        messge = ' $'
*        goto 999
*      endif
**
** set flags to indicate stacking type
*      recrsv    = .false.
*      xplcit    = .false.
*      rndm      = .false.
*      inf_thick = .false.
** initialize l_cnt
*      l_cnt = 0
*      if(iflag.eq.1) then
*        xplcit = .true.
*      else if(iflag.eq.2) then
*        recrsv = .true.
*      endif
**
*      if(recrsv) then
*        messge = 'Problem reading crystal thickness for recursion.$'
*        call GETLNE(unit_no, line, *999)
*        call TOUPPR(line)
**
*        if(index(line,'INFINITE').eq.0) then
*          inf_thick = .false.
*          messge =
*     |       'Problems reading the number of RECURSIVE layers.$'
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) line
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) l_cnt
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**          read(line, *, err=999) l_cnt
** check l_cnt
*          if(l_cnt.eq.0) then
*            messge = 'The number of RECURSIVE layers cannot be zero.$'
*          else if(l_cnt.lt.0) then
*            messge =
*     |           'The number of RECURSIVE layers cannot be negative.$'
*          else if(l_cnt.gt.RCSV_MAX) then
*            write(op,104)
*     |       'WARNING: Number of RECURSIVE layers ', l_cnt,
*     |       ' exceeds the maximum ', RCSV_MAX
*            write(op,105) 'An INFINITE thickness is assumed.'
** re-set l_cnt and proceed assuming infinite thickness
*            l_cnt = 0
*            inf_thick = .true.
*          endif
*        else
*          inf_thick = .true.
*        endif
**
*      else if(xplcit) then
** read in layer sequence
*        messge = 'Problem reading EXPLICIT layer sequencing.$'
*        call GETLNE(unit_no, line, *999)
*        call TOUPPR(line)
*        rndm = index(line,'RANDOM ').ne.0
** catch possible error - user requested INFINITE
*        if(.not.rndm) then
*          if(index(line,'INFINITE').ne.0) then
*            write(op,107) 'Maximum number of layers allowed is ',XP_MAX
*            messge = 'Too many layers for an EXPLICIT sequence.$'
*            goto 999
*          endif
*        endif
**
*        if(rndm) then
** the user asked for a random layer sequencing.
**
*          messge =
*     |       'Problems reading the number of random EXPLICIT layers.$'
*          write(tmpline, '(a)', err=999)
*     |        line(LENGTH('RANDOM')+2:len(line))
** for strict FORTRAN77 compliance
*          write(scrtch, *, err=990) tmpline
*          rewind(scrtch, err=990)
*          read(scrtch, *, err=999) l_cnt
*          rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**          read(tmpline, *, err=999) l_cnt
**
*        else
*          l_cnt = 1
*          messge =
*     |       'Problems reading EXPLICIT layer sequence in RDSTAK.$'
*   20     i = LAYCNT(line)
**
*          if(i.ge.0 .and. i.le.40) then
**
** for strict FORTRAN77 compliance
*            write(scrtch, *, err=990) line
*            rewind(scrtch, err=990)
*            read(scrtch, *, err=999) (l_seq(j), j=l_cnt,l_cnt+i-1)
*            rewind(scrtch, err=990)
** for more lenient, efficient, compilers
**            read(line, *, err=999) (l_seq(j), j=l_cnt,l_cnt+i-1)
**
*            l_cnt = l_cnt + i
*          else
*            goto 999
*          endif
**
*          call GETLNE(unit_no, line, *999)
*          call TOUPPR(line)
** check we have not run into the 'TRANSITIONS' section yet.
*          if(index(line,'TRANSITIONS').ne.1) goto 20
*          l_cnt = l_cnt - 1
**
** If we are here, we have gone one line too far. Back up.
*          write(messge,103)
*     |       'Problems ''backspacing'' unit ', unit_no, '.$'
*          backspace(unit = unit_no, err = 980)
**
*        endif
**
** check we didn't enter too many layers
*        if(l_cnt.gt.XP_MAX) then
*          write(op,104) 'WARNING: No. of EXPLICIT layers, ',l_cnt,
*     |      ' exceeds maximum of ',XP_MAX
*          write(op,106) XP_MAX, ' EXPLICIT layers assumed'
*          l_cnt = XP_MAX
*        endif
**
** check that all the layer numbers are within bounds
*        if(.not.rndm) then
*          do 30 i = 1, l_cnt
*            if(l_seq(i).gt.n_layers) then
*              write(messge,102) 'Illegal layer number,', l_seq(i),
*     |         '. Number must not exceed', n_layers, '.$'
*              goto 980
*            else if(l_seq(i).lt.1) then
*              write(messge,103) 'Illegal layer number,', l_seq(i),
*     |         '. Number must be greater than 0.$'
*              goto 980
*            endif
*   30     continue
*        endif
*      endif
**
*      RDSTAK = .true.
*  980 return
*  990 messge = 'Problems using scratch file in RDSTAK.$'
*  999 call WRTLNE(line)
*      return
*  100 format(1x, 'ERROR: ', a)
*  101 format(1x, 2a)
*  102 format(2(a, i4), a)
*  103 format(a, i2, a)
*  104 format(1x, 2(a, i4))
*  105 format(1x, a)
*  106 format(1x, i4, a)
*  107 format(1x, a, i5)
*      end
**
** ______________________________________________________________________
** Title: RDTRNS
** Authors: MMJT
** Date: 4 Oct 1989; 15 Mar 1995; 29th Aug 2000
** Description: This function reads the layer stacking probabilities
** and stacking vectors.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**            uses:  n_layers, recrsv, rndm, xplcit
**
**        modifies:  l_alpha, l_r, r_B11, r_B22, r_B33, r_B12, r_B23,
**                   r_B31
**
**      RDTRNS returns logical .true. if all the stacking
**      parameters were properly specified.
** ______________________________________________________________________
**
*      logical function RDTRNS(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line, tmpline, arg
*      integer*4 PARENT, NXTARG, i, j, n, tempfact
*      real*8 sum, RDNMBR
*      logical*4 ok
**
** external function
*      external PARENT, NXTARG, RDNMBR
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR, WRTLNE
**
*      RDTRNS = .false.
**
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
*      if(index(line,'TRANSITIONS ').ne.1) then
*        messge = '''TRANSITIONS'' section not found.$'
*        goto 999
*      endif
**
*      do 10 i = 1, n_layers
*        do 20 j = 1, n_layers
*          write(messge,100)
*     |     'Unexpected EOF at vector ',i,j,' in ''TRANSITIONS''.$'
*          call GETLNE(unit_no, line, *999)
** Track down if there are parentheses in this line. If so then there
** are FatsWaller (Debye-Waller for layers!) factors
*          write(messge,100)
*     |      'Confusing use of parentheses in data for vector ',i,j,'.$'
*          tempfact = PARENT(line)
*          if(tempfact.lt.0) goto 999
** initialize
*          l_r(1,j,i) = ZERO
*          l_r(2,j,i) = ZERO
*          l_r(3,j,i) = ZERO
*
*          r_B11(j,i) = ZERO
*          r_B22(j,i) = ZERO
*          r_B33(j,i) = ZERO
*          r_B12(j,i) = ZERO
*          r_B31(j,i) = ZERO
*          r_B23(j,i) = ZERO
**
*          tmpline = line
**
*          if(tempfact.eq.1) then
*            write(messge,100)
*     |     'Expected alpha_ij, R_ij, (dR_ij) for vector ',i,j,'.$'
*          else
*            write(messge,100)
*     |     'Expected alpha_ij, R_ij data for vector ',i,j,'.$'
*          endif
**
** read stacking probability
*          n = NXTARG(tmpline,arg)
*          if(n.le.0) goto 999
*          l_alpha(j,i) = RDNMBR(arg,ok)
*          if(.not.ok) goto 999
** If zero, there is no need to read further
*          if(dabs(l_alpha(j,i)).lt.eps6) goto 20
**           
** read stacking x-vector
*          n = NXTARG(tmpline,arg)
*          if(n.le.0) goto 999
*          l_r(1,j,i) = RDNMBR(arg,ok)
*          if(.not.ok) goto 999
**           
** read stacking y-vector
*          n = NXTARG(tmpline,arg)
*          if(n.le.0) goto 999
*          l_r(2,j,i) = RDNMBR(arg,ok)
*          if(.not.ok) goto 999
**           
** read stacking z-vector
*          n = NXTARG(tmpline,arg)
*          if(n.le.0) goto 999
*          l_r(3,j,i) = RDNMBR(arg,ok)
*          if(.not.ok) goto 999
**
** Read temperature factors, if present
*          if(tempfact.eq.1) then
**           
** read B11 "Fats-Waller"
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            r_B11(j,i) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**           
** read B22 "Fats-Waller"
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            r_B22(j,i) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**           
** read B33 "Fats-Waller"
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            r_B33(j,i) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**           
** read B12 "Fats-Waller"
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            r_B12(j,i) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**           
** read B31 "Fats-Waller"
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            r_B31(j,i) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**           
** read B23 "Fats-Waller"
*            n = NXTARG(tmpline,arg)
*            if(n.le.0) goto 999
*            r_B23(j,i) = RDNMBR(arg,ok)
*            if(.not.ok) goto 999
**           
*          endif
**
** stacking probabilities should not be negative.
*          if(l_alpha(j,i).lt.ZERO) then
*            write(messge,100) 'alpha',i,j,' is negative.$'
*            goto 999
*          endif
** the Rz(i,i) vectors should not be zero or negative.
** if an i-i transition exists
*          if(j.eq.i) then
*            if(l_alpha(i,i).ne.ZERO) then
*              if(l_r(3,i,i).le.ZERO) then
*                write(messge,100)
*     |'Invalid negative (or zero) z-component for Rz',i,i,'.$'
*                goto 999
*              endif
*            endif
*          endif
*   20   continue
*   10 continue
**
** Check stacking probabilities.
** Generate error message if they do not sum to 1, and we are using
** a RECURSIVE stacking description. If EXPLICIT stacking, then we
** do not use the probabilities (but we still assume they are there
** to be read) but generate a warning.
**
*      do 30 i = 1, n_layers
*        sum = ZERO
*        do 40 j = 1, n_layers
*          sum = sum + l_alpha(j,i)
*   40   continue
*        if(abs(sum - ONE).gt.eps6) then
*          write(messge,101)
*     |     'Stacking probabilities from LAYER ',i,' do not sum to 1.$'
*          if(recrsv.or.rndm) then
*            goto 980
*          else if(xplcit .and. .not.rndm) then
*            write(op,102) messge(1:index(messge, '$') - 1)
*          endif
*        endif
*   30 continue
**
*      RDTRNS = .true.
*  980 return
*  999 call WRTLNE(line)
*      return
*  100 format(a, '(', i2, ',', i2, ')', a)
*  101 format(a, i2, a)
*  102 format(1x, 'WARNING: ', a)
*      end
**
** ______________________________________________________________________
** Title: RDWAVL
** Authors: MMJT
** Date: 4 Oct 1989
** Description: This function reads the radiation wavelength.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**        modifies:  lambda
**
**      RDWAVL returns logical .true. if an acceptable value
**      for the radiation wavelength was read.
** ______________________________________________________________________
**
*      logical function RDWAVL(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line
**
** external subroutines (Some compilers need them declared external)
**      external GETLNE, WRTLNE
**
*      RDWAVL = .false.
**
*      messge = 'lambda is specified improperly.$'
*      call GETLNE(unit_no, line, *999)
**
** for strict FORTRAN77 compliance
*      write(scrtch, *, err=990) line
*      rewind(scrtch, err=990)
*      read(scrtch, *, err=999) lambda
*      rewind(scrtch, err=990)
**
** for more lenient, efficient compilers
**      read(line, *, err=999) lambda
**
*      if(lambda.le.ZERO) then
*        messge = 'lambda is negative (or zero).$'
*        goto 999
*      endif
**
*      RDWAVL = .true.
*      return
*  990 messge = 'Problems using scratch file in RDWAVL.$'
*  999 call WRTLNE(line)
*      return
*      end
**
** ______________________________________________________________________
** Title: RDWDTH
** Authors: MMJT
** Date: 3 Mar 1995
** Description: This function reads the lateral layer dimensions 
** from the data file. First we check that the data on the line does 
** not refer to the layer atomic positions. If it does, assume infinite
** widths, backspace the file, and return to let RDLAYR do its job.
** Unlike the other read functions, RDWDTH does not return an error in
** this instance, so as to maintain compatibility with data files 
** written for v 1.7xx which did not treat lateral layer width
** broadening.
**
**      ARGUMENTS:
**            unit_no  -  Logical unit number that the data file is to
**                        be read from. (input).
**            messge   -  A short message indicating what went wrong
**                        (if anything) during the datafile read.
**                        messge is terminated by a token '$'. (output).
**
**      COMMON VARIABLES:
**        modifies:  finite_width, Wa, Wb
**
**      RDWAVL returns logical .true. if no errors were encountered.
** ______________________________________________________________________
**
*      logical function RDWDTH(unit_no, messge)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      integer*4 unit_no
*      character*(*) messge
**
*      character*200 line
**
*      integer*4 n, CNTARG
**
** external function
*      external CNTARG
** external subroutines (Some compilers need them declared external)
**      external GETLNE, TOUPPR, WRTLNE
**
*      RDWDTH = .false.
**
** Does the next line belong to the layer structural data? If so, 
** assume infinite layer widths. No error message is generated so 
** that there is compatibility with data files for versions 1.7xx.
*      messge = 'layer width, or layer data, are specified improperly.$'
*      call GETLNE(unit_no, line, *999)
*      call TOUPPR(line)
*      if(index(line, 'LAYER').eq.1) then
*        finite_width = .false.
*        write(messge,100)
*     |   'Problems ''backspacing'' unit ',unit_no, '.$'
** backspace file so that RDLAYR will work properly
*        backspace(unit = unit_no, err = 990)
*        goto 900
*      endif
**
** assume that layer width data is specified
*      if(index(line,'INFINITE').eq.1) then
*        finite_width = .false.
*        goto 900
*      else
*        finite_width = .true.
*      endif
**
** are both Wa and Wb specified, or just a generic width?
*      messge = 'layer characteristic widths are specified improperly.$'
*      n = CNTARG(line)
**
** for strict FORTRAN77 compliance
*      write(scrtch, *, err=990) line
*      rewind(scrtch, err=990)
**
*      if(n.eq.1) then
** for strict FORTRAN77 compliance
*        read(scrtch, *, err=999) Wa
*        rewind(scrtch, err=990)
** for more lenient, efficient compilers
**        read(line, *, err=999) Wa
*        Wb = Wa
*      else if(n.eq.2) then
** for strict FORTRAN77 compliance
*        read(scrtch, *, err=999) Wa, Wb
*        rewind(scrtch, err=990)
** for more lenient, efficient compilers
**        read(line, *, err=999) Wa, Wb
*      else
*        goto 999
*      endif
**
*      if(finite_width) then
*        if(Wa.gt.inf_width .and. Wb.gt.inf_width) then
*          finite_width = .false.
*          write(op,101) 'Layers will be treated as if infinitely wide.'
*        endif
*        if(Wa.le.ZERO .or. Wb.le.ZERO) then
*          write(op,101)'Illegal values for layer characteristic widths'
*          goto 999
*        endif
*      endif
**
*  900 RDWDTH = .true.
*      return
*  990 messge = 'Problems using scratch file in RDWDTH.$'
*  999 call WRTLNE(line)
*      return
*  100 format(a, i2, a)
*  101 format(1x, a)
*      end
**
** ______________________________________________________________________
* Title: RAN3
* Authors: Press, Flannery, Teukolsky and Vetterling
* Date: Copyright (C) 1985
* Returns a uniform random deviate between 0.0 and 1.0. Set 'idum'
* to any negative value to initialize or reinitialize the sequence.
* This version is modified to return real*8 values, and enforces static
* storage of all local variables by use of the 'save' statement
* (In fact 'seed' is the important variable to save, but we save all
* anyway).
*
*      ARGUMENTS:
*            idum       -  Set -ve to initialize. (input).
*
*      RAN3 returns a real random number between 0 and 1
* ______________________________________________________________________
*
      real*8 function RAN3(idum)
      implicit none
      save
*
      integer*4 idum
*
      real*8 big, seed, mz, fac
      parameter (big=4000000.0D0,seed=1618033.0D0,mz=0.0D0,fac=2.5D-7)
      real*8 ma(55)
      real*8 mj, mk
      integer*4 iff, ii, i, j, inext, inextp
*
      data iff /0/
*
      if(idum.lt.0 .or. iff.eq.0) then
        iff = 1
        mj = seed - iabs(idum)
        mj = mod(mj,big)
        ma(55) = mj
        mk = 1.0D0
        do 11 i = 1, 54
          ii = mod(21*i,55)
          ma(ii) = mk
          mk = mj - mk
          if(mk.lt.mz) mk = mk + big
          mj = ma(ii)
   11   continue
        do 13 j = 1, 4
          do 12 i = 1, 55
            ma(i) = ma(i) - ma(1 + mod(i + 30,55))
            if(ma(i).lt.mz) ma(i) = ma(i) + big
   12     continue
   13   continue
        inext = 0
        inextp = 31
        idum = 1
      endif
*
      inext = inext + 1
      if(inext.eq.56) inext = 1
      inextp = inextp + 1
      if(inextp.eq.56) inextp = 1
      mj = ma(inext) - ma(inextp)
      if(mj.lt.mz) mj = mj + big
      ma(inext) = mj
      RAN3 = mj * fac
*
      return
      end
*
* ______________________________________________________________________
* Title: SALUTE
* Author: MMJT
* Date: 19th May, 2010
* Draws the DIFFaX flag. Saluting is optional.
*
*      ARGUMENTS:
*           No arguments are needed.
* ______________________________________________________________________
*
      subroutine SALUTE()
      include 'DIFFaX.par'
*     save
*
      write(op,1) '***************************************************'
      write(op,1) '***************************************************'
      write(op,1) '*                                                 *'
      write(op,1) '*  DDDD     II   FFFFFF   FFFFFF           X   X  *'
      write(op,1) '*  D   D    II   F        F        aaaa     X X   *'
      write(op,1) '*  D    D   II   FFFF     FFFF    a    a     X    *'
      write(op,1) '*  D   D    II   F        F       a   aa    X X   *'
      write(op,1) '*  DDDD     II   F        F        aaa a   X   X  *'
      write(op,1) '*                                                 *'
      write(op,1) '***************************************************'
      write(op,1) '****************** DIFFaX v1.813 ******************'
      write(op,1) '***************************************************'
      write(op,1) '***************** 19th May,  2010 *****************'
      write(op,1) '***************************************************'
      write(op,1) '*                                                 *'
      write(op,1) '*   A computer program for calculating            *'
      write(op,1) '*   Diffraction Intensity From Faulted Crystals   *'
      write(op,1) '*                                                 *'
      write(op,1) '*   Authors: Michael M. J. Treacy                 *'
      write(op,1) '*            Michael W. Deem                      *'
      write(op,1) '*                                                 *'
      write(op,1) '***************************************************'
      write(op,1)
*
      return
    1 format(1x, a51)
      end
*
* ______________________________________________________________________
** Title: SFC
** Authors: MWD and MMJT
** Date: 11 Feb 90; 7 Mar 1995
** Description: This routine reads the atomic scattering factor constants
** from the file 'sfname' (usually given as "data.sfc"). The scattering
** factors are used to generate the layer form factors. SFC returns
** .false. if there are too many distinct atom types (ie. more than
** MAX_TA). The value of MAX_TA can be reset in 'DIFFaX.par'.
**
**      ARGUMENTS:
**           No arguments are used. All data is in 'COMMON'.
**
**      COMMON VARIABLES:
**            uses:  n_actual, l_n_atoms, a_name, n_atoms, rad_type
**                   NEUTRN, X_RAY, rad_type, sfname
**
**        modifies:  n_atoms, x_sf, n_sf, e_sf, a_type
**
**      SFC returns logical .true. if it read the scattering factor data
**      safely.
** ______________________________________________________________________
**
*      logical function SFC()
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      logical ok, new_atom, our_atom, done
*      logical list(MAX_TA+1)
*      integer*4 i, j, n, m, LENGTH
*      character*4 name
*      character*120 line
*      real*4 tmp
**
** external subroutine (Some compilers need them declared external)
**      external ATOMS
** external function
*      external LENGTH
**
*      write(op,301) 'Reading scattering factor datafile ''',
*     |                   sfname(1:LENGTH(sfname)),'''. . .'
**
*      SFC = .false.
*      do 10 i = 1, MAX_TA
*        list(i) = .false.
*        atom_l(i) = ' '
*   10 continue
*      m = 0
** determine all distinct atom types
*      do 30 i = 1, n_actual
*        do 40 j = 1, l_n_atoms(i)
*          name = a_name(j,i)
*          call TOUPPR(name)
** see if this name has been seen already
*          n = 0
*   50     n = n + 1
*            new_atom = name.ne.atom_l(n)
*            if(new_atom .and. n.lt.m) goto 50
*          if(new_atom) then
*            m = m + 1
*            if(m.gt.MAX_TA) goto 220
*            atom_l(m) = name
*            a_type(j,i) = m
*          else
*            a_type(j,i) = n
*          endif
*   40   continue
*   30 continue
*      n_atoms = m
** now find data for each atom type in file
** pass through file only once
*   60 read(sf, 300, end=90, err=200) line
*        name = line(1:4)
*        call TOUPPR(name)
** see if this is one of the distinct atoms
*        i = 0
*   70   i = i + 1
*          our_atom = name.eq.atom_l(i)
*        if(.not.our_atom .and. i.lt.n_atoms) goto 70
** have we read all that we need to know?
*        done = .true.
*        do 80 n = 1, n_atoms
*          done = done.and.list(n)
*   80   continue
** If we're done, close file.
*        if(done) goto 90
*        if(our_atom .and. .not.list(i)) then
** mark this atom's data as read in
*          list(i) = .true.
** and read it in
*          if(rad_type.eq.X_RAY) then
*            read(line,310,err=210) (x_sf(j,i), j=1, 9)
*          else if(rad_type.eq.NEUTRN) then
*            read(line,320,err=210) n_sf(i)
*          else
*            read(line,310,err=210) (x_sf(j,i), j=1, 9), tmp, e_sf(i)
*          endif
*        endif
*      goto 60
*   90 close(sf,err=240)
** see if all the data for each atom has been read in
*      ok = .true.
*      do 100 i = 1, n_atoms
*        if(.not.list(i)) then
*          ok = .false.
*          write(op,330) 'ERROR: Data for atom ''', atom_l(i),
*     |     ''' not found in file ''', sfname(1:LENGTH(sfname)), ''''
*          call ATOMS()
*        endif
*  100 continue
*      if(ok) then
*        SFC = .true.
*        write(op,400) 'Scattering factor data read in.'
*      endif
*      return
*  200 write(op,301) 'Scattering factor file ''',
*     |    sfname(1:LENGTH(sfname)), ''' defective.'
*      close(sf,err=240)
*      return
*  210 write(op,400) 'ERROR reading scattering factor data.'
*      close(sf,err=240)
*      return
*  220 write(op,402) 'There are too many types of atoms in layer ', i
*      write(op,400) 'Atoms recorded so far are:'
*      do 110 j = 1, MAX_TA
*        write(op,403) '       Atom type ', j, '   ', atom_l(j)
*  110 continue
*      write(op,402) '  Maximum number of types allowed is ', MAX_TA
*      return
*  240 write(op,301) 'Unable to close scattering factor file ''',
*     |    sfname(1:LENGTH(sfname)), '''.'
*      return
*  300 format(a)
*  301 format(1x, 3a)
*  310 format(t5, 10f11.6, i3)
*  320 format(t104, f11.6)
*  330 format(1x, 5a)
*  400 format(1x, a)
*  401 format(1x, 2a)
*  402 format(1x, a, i2)
*  403 format(1x, a, i2, 2a)
*      end
**
** ______________________________________________________________________
* Title: SHARP
* Author: MMJT
* Date: 29 Aug 1991
* Description: This subroutine determines whether or not
* spots on a given h, k row are sharp. It does this by examining
* the intensity at l = 0 (or, if absent, at l = d_l) and comparing
* with the intensity at a nearby l-value. If the peak is sharp, there
* will be a large change in intensity.
*
*      ARGUMENTS:
*            h    -  reciprocal vector h-component. (input).
*            k    -  reciprocal vector k-component. (input).
*            d_l  -  steps in l-component of the reciprocal lattice
*                    vector that sharp spots are likely to be found at.
*                                                            (input).
*
*      COMMON VARIABLES:
*            uses:  a0, b0, c0, d0, lambda, PI
*
*        modifies:  no COMMON variables are modified
*
*      SHARP returns logical .true. if it thinks h, k contains sharp
*      spots.
* ______________________________________________________________________
*
      logical function SHARP(h, k, d_l)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k
      real*8 d_l
*
      logical ok
      real*8 S, PNTINT, ANGLE, LL, l, i1, i2, x, theta, l_next
*
* external subroutine (Some compilers need them declared external)
*      external GET_F
* external functions
      external PNTINT
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
* ANGLE is the Bragg angle (in radians) of the h,k,l plane
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
* LL is the maximum allowable l value for a given h,k and theta
      LL(theta,h,k) =  sqrt( ( (TWO*sin(theta)/lambda)**2
     |                         - h*h*a0 - k*k*b0 - h*k*d0) / c0 )
*
      SHARP = .false.
*
* get the intensity at hkl, with l = 0 initially
      l = ZERO
   10 i1 = PNTINT(h, k, l, ok)
      if(.not.ok) goto 100
* If there is an extinction at hkl, try again at l = l + d_l
      if(i1.lt.eps4) then
        l = l + d_l
        if(ANGLE(h,k,l).lt.HALF*PI) then
          goto 10
        else
          goto 100
        endif
      endif
*
* Define a spot to be sharp if intensity is halved at l = l + d_l/100
      theta = ANGLE(h, k, l)
      x = min(d_theta, HALF*th2_max-theta)
      l_next = LL(theta+x, h, k)
      i2 = PNTINT(h, k, l+l_next*eps2, ok)
      if(.not.ok) goto 100
*
      SHARP = i1 .gt. TWO*i2
*
  100 return
      end
*
* ______________________________________________________________________
* Title: SMUDGE
* Author: MMJT
* Date: 30 Oct 1989; 16 April 1999
* Description: This subroutine convolutes 'array', of length
* 'arrsize' with a Gaussian of standard deviation 'sigma'.
* NOTE: This routine does not conserve area under the convoluted curve.
* It is designed so that the peak heights are unchanged.
*
*      ARGUMENTS:
*            array    -  The name of the input array. (input).
*            arrsize  -  Size of the input array. (input).
*            sigma    -  Standard deviation of the Gaussian as a
*                        multiple of the array sampling step. (input).
*            ok       -  logical flag indicating all went well.
*                                                      (output).
* ______________________________________________________________________
*
      subroutine SMUDGE(array, arrsize, sigma, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*     save
*
      integer*4 arrsize
      real*8 array(arrsize), sigma
      logical ok
*
      integer*4 m, i, j
      real*8 tmparr(SADSIZE)
      real*8 k1, k2, tmp, tmp1, tmp2, gss, normalize
*
      if(sigma.eq.ZERO) return
*
      do 10 i = 1, arrsize
        if(array(i).gt.maxsad) array(i) = maxsad
        tmparr(i) = array(i)
        array(i) = ZERO
   10 continue
*
* go out to 5 standard deviations, or to the end of the spectrum
      m = nint(FIVE * sigma)
      k1 = HALF / ( sigma * sigma )
      if(m.gt.arrsize) m = arrsize
* get normalization constant. We wish peak heights to be unchanged.
      normalize = ONE
      do 20 i = 1, m
        normalize = normalize + TWO * exp( -k1 * dble(i*i) )
   20 continue
*
      if(normalize.eq.ZERO) then
        write(op,100) 'ERROR in SMUDGE: Zero normalization constant.'
        ok = .false.
        goto 999
      endif
      normalize = ONE / normalize
*
      do 30 i = 0, m
        k2 = k1*dble(i*i)
        gss = exp(-k2)
        do 40 j = 1, arrsize
          tmp1 = ZERO
          tmp2 = ZERO
          if((j-i).gt.0) tmp1 = tmparr(j-i)
          if((j+i).le.arrsize) tmp2 = tmparr(j+i)
          tmp = tmp1 + tmp2
          if(i.eq.0) tmp = HALF * tmp
          array(j) = array(j) + gss * tmp * normalize
   40   continue
   30 continue
*
  999 return
  100 format(1x, a)
      end
*
* ______________________________________________________________________
* Title: SPHCST
* Author: MWD
* Date: 18 Aug 1988
* Description: This subroutine determines the constants used in
* determining the magnitude of reciprocal lattice vectors.
*
*      ARGUMENTS:
*            No input arguments.
*
*      COMMON VARIABLES:
*            uses:  cell_a, cell_b, cell_c, cell_gamma
*
*        modifies:  a0, b0, c0, d0, ab0, bc0, ca0
* ______________________________________________________________________
*
      subroutine SPHCST()
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      a0 = ONE / (cell_a * sin(cell_gamma))**2
      b0 = ONE / (cell_b * sin(cell_gamma))**2
      c0 = ONE / (cell_c)**2
      d0 = -TWO * cos(cell_gamma) /
     |            (cell_a * cell_b * sin(cell_gamma)**2)
*
      ab0 = sqrt(a0 * b0)
      bc0 = sqrt(b0 * c0)
      ca0 = sqrt(c0 * a0)
*
      return
      end
*
* ______________________________________________________________________
** Title: STREAK
** Author: MWD & MMJT
** Date: 18 Aug 1988, revised 22 June 1989.
** Description:  This routine outputs the integrated (ie. averaged)
** intensitites from h,k,l0 to h,k,l1, in steps of dl, to the file
** strkfile.
**
**      ARGUMENTS:
**            FN        -  Function name passed by reference. The
**                         choice is between GLQ16 (non-adaptive
**                         Gauss-Legendre integration), and AGLQ16
**                         (adaptive Gauss-Legendre integration).
**                                                            (input).
**            strkfile  -  The name of the file that streak data is
**                         to be output to. (input).
**            ok        -  logical flag indicating all went well.
**                                                      (output).
**
**      COMMON VARIABLES:
**            uses:  a0, b0, c0, d0, lambda, cntrl, CFile, xplcit,
**                   rad_type, X_RAY
**
**        modifies:  no COMMON variables are modified
** ______________________________________________________________________
**
*      subroutine STREAK(FN, strkfile, ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      real*8 FN
*      character*(*) strkfile
*      logical ok
**
*      logical its_hot
*      integer*4 h, k, LENGTH, i, i_step
*      real*8 l, theta, x, ANGLE, S, Q2, l0, l1
*      real*8 dl, W4, intervals
*      parameter (intervals = TWENTY)
**
** external functions
*      external FN, LENGTH
** external subroutines (Some compilers need them declared external)
**      external XYPHSE, PRE_MAT, GET_F
**
** statement functions
** S is the value of 1/d**2 at hkl
*      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
** ANGLE is the Bragg angle (in radians) of the h,k,l plane
*      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
** W4 is the X-ray polarization factor
*      W4(theta) = HALF * (ONE + (cos(TWO*theta))**2)
**
*      Q2 = FOUR / (lambda**2)
*   10 write(op,400) 'Enter h, k, l0, l1, delta l'
*      read(cntrl,*,err=10) h, k, l0, l1, dl
*      if(CFile) write(op,401) h, k, l0, l1, dl
** check input
*      if(l1.eq.l0) then
*        write(op,400) 'Illegal input: l0 equals l1'
*        goto 999
*      else if(dl.eq.ZERO) then
*        write(op,400) 'Illegal zero value of dl entered'
*        write(op,402)'A value of ',(l1-l0)/(FIVE*HUNDRED),' is assumed'
*      else if(l1.gt.l0 .and. dl.lt.ZERO) then
*        write(op,400) 'l1 is greater than l0. +ve dl assumed'
*        dl = -dl
*      else if(l1.lt.l0 .and. dl.gt.ZERO) then
*        write(op,400) 'l0 is greater than l1. -ve dl assumed'
*        dl = -dl
*      endif
** The origin may be hotter than hell! Let's check first.
*      its_hot = (h.eq.0 .and. k.eq.0) .and.
*     |         l0*l1.le.ZERO .and. rad_type.eq.ELECTN
*      if(its_hot) then
*        write(op,400) 'Cannot scan the origin with electron radiation'
*        write(op,400) 'Origin will be skipped.'
*      endif
** check angles are meaningful
*      if(S(h,k,l0).gt.Q2 .or. S(h,k,l1).gt.Q2) then
*        if(S(h,k,l0).gt.Q2) write(op,403) h, k, l0,
*     |            ' exceeds 180 degree scattering angle!'
*        if(S(h,k,l1).gt.Q2) write(op,403) h, k, l1,
*     |            ' exceeds 180 degree scattering angle!'
*        goto 10
*      endif
**
*      write(op,404) 'Writing streak data to file ''',
*     |      strkfile(1:LENGTH(strkfile)),'''. . .'
*      call XYPHSE(h, k)
*      call PRE_MAT(h, k)
**
*      i_step = nint( (l1 - l0) / (dl * intervals) )
** Immune system to the rescue!
*      if(i_step.le.0) i_step = 10
*      i = 0
*      do 20 l = l0, l1, dl
** If we are dealing with electrons, make sure we avoid the origin
*        if(its_hot .and. l*(l+dl).le.ZERO) then
*          x = ZERO
*          goto 30
*        endif
*        i = i + 1
*        if(mod(i,i_step).eq.0) write(op,405) 'Reached l = ',l
*        x = FN(h,k,l,l+dl,ok)
*        if(.not.ok) goto 999
** note: since this is streak data, only the X_RAY input needs
** correcting for polarization.
*        if(rad_type.eq.X_RAY)  x = x * W4(ANGLE(h,k,l+HALF*dl))
*   30   write(sk,406,err=100) l, char(9), x
*   20 continue
*      if(sk.ne.op) close(sk,err=110)
*      write(op,404) 'Streak data file, ''',
*     |      strkfile(1:LENGTH(strkfile)),''' written to disk.'
*      return
*  100 write(op,404) 'ERROR writing to streak data file ''',
*     |      strkfile(1:LENGTH(strkfile)),''''
*      if(sk.ne.op) close(sk,err=110)
*      return
*  110 write(op,404) 'Unable to close streak data file ''',
*     |      strkfile(1:LENGTH(strkfile)),''''
*      return
*  999 write(op,405) 'ERROR encountered in streak integration at l = ',l
*      return
*  400 format(1x, a)
*  401 format(1x, 2i3, 3f10.5)
*  402 format(1x, a, f10.5, a)
*  403 format(1x, 2i3, f10.5, a)
*  404 format(1x, 3a)
*  405 format(1x, a, f10.5)
*  406 format(1x, e12.5, a, e14.6)
*      end
**
** ______________________________________________________________________
* Title: THRESH
* Author: MMJT
* Date: 21 Jan 1995
* Samples intensities in reciprocal space to get a "feeling" for what
* kind of values are out there prior to testing diffraction symmetry.
* CHK_SYM and GET_SYM measure the fractional deviations of intensity
* from (potentially) symmetry-related points in reciprocal space.
* This method runs into problems when the intensity being measured
* is close to zero. Miniscule intensity variations can appear to be
* huge relative to zero!
* This function is needed in order to obtain a (crude) estimate of
* which intensity values are too small to worry about, even if the 
* relative intensity variations seem to be large.
* This function will be of no use if there are no streaks, that is,
* if the crystal is perfect.
*
*      ARGUMENTS:
*            ok       -  logical flag indicating all went well.
*                                                      (output).
*
*      COMMON VARIABLES:
*            uses:  a0, b0, c0, d0, lambda, no_trials,
*                   max_angle
*
*        modifies:  ok, max_angle, h_bnd, k_bnd, l_bnd, tiny_inty
*
* ______________________________________________________________________
*
      subroutine THRESH(ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      logical ok
*
      integer*4 i, h, k, idum
      real*8 RAN3, PNTINT, S, ANGLE
      real*8 l, tot_int
*
* external functions
      external RAN3, PNTINT
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
*
* initialize random numbers in RAN3
      idum = -1
*
* First define angular range to sample. h_bnd, k_bnd, l_bnd
* (defined in HKL_LIM) and max_angle, are used later on 
* in GET_SYM and CHK_SYM.
      max_angle = QUARTER * PI
      call HKL_LIM()
*
* Sample the typical intensities that are out there
      tot_int = ZERO
      do 10 i = 1, no_trials
   20   h = int( dble(h_bnd + 1) * RAN3(idum) )
        k = int( dble(k_bnd + 1) * RAN3(idum) )
        l = l_bnd * RAN3(idum)
* make sure we are not sampling at too high an angle
        if(TWO * ANGLE(h, k, l) .gt. max_angle) goto 20
* get I(h,k,l)
        tot_int = tot_int + PNTINT(h, k, l, ok)
        if(.not.ok) goto 999
   10 continue
*
* Estimate some suitable fraction of the average intensity. This
* fraction defines a baseline intensity equivalent to zero for
* the purposes of the symmetry testing later on.
*
      tiny_inty = tot_int * eps5 / no_trials
*
  900 return
  999 write(op,100) 'ERROR in intensity calculation in THRESH'
      write(op,200) '   at h,k,l = ', h,',', k,',', l
      return
  100 format(1x, a)
  200 format(1x, a, i3, a, i3, a, f7.2)
      end
*
* ______________________________________________________________________
** Title: TOUPPR
** Author: MMJT
** Date: 3 August 1989
** Description: Converts alphabetic characters in 'line' to uppercase.
** Works for ascii characters, since it assumes alphabetic characters
** are contiguous.
** The neat trick of using
**                 lower_case = ichar(chr).and.95
** will not always work in ansi77 fortran, as some compilers evaluate
** the right hand side as a logical .true. or .false. rather than a
** bitwise 'anded' integer.
** The alternative
**                 lower_case = ichar(chr).iand.95
** is equally undesirable since the operator .iand. is a military
** extension, and is not ansi77.
**
**      ARGUMENTS:
**            line  -  Input line of characters to be converted to
**                     upper case. The converted string is returned
**                     in 'line'.
** ______________________________________________________________________
**
*      subroutine TOUPPR(line)
**
*      character*(*) line
**
*      integer*4 a, z, offset, i, itmp, lin_len
**
*      a = ichar('a')
*      z = ichar('z')
*      offset = a - ichar('A')
*      lin_len = len(line)
**
*      i = 1
*   10 itmp = ichar(line(i:i))
*        if(itmp.ge.a .and. itmp.le.z) line(i:i) = char(itmp - offset)
*        i = i + 1
*        if(i.le.lin_len)
*     |goto 10
**
*      return
*      end
**
** ______________________________________________________________________
* Title: TRMSPC
* Author: MMJT
* Date: 20 April 1989; 7 Mar 1995
* Description:  This function locates a suitable cut-off angle,
* below which we can ignore the huge intensities which may occur
* close to the origin when the full adaptive integration option is used.
* TRMSPC is called only when theta = 0 is in the integrated range.
*
*      ARGUMENTS:
*            th2_low  -  a 2theta cut-off angle to be found (output)
*
*      COMMON VARIABLES:
*            uses:  th2_min, th2_max, d_theta, spec, RAD2DEG
*
*        modifies:  spec
*
*            TRMSPC returns logical .true. if all went well.
* ______________________________________________________________________
*
      logical function TRMSPC(th2_low)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      real*8 th2_low
*
      integer*4 i, i_min, i_max
*
      i_max = int(HALF*(th2_max - th2_min) / d_theta) + 1
* spec(1) corresponds to the intensity at the origin and is always zero.
      i = 2
   20 i = i + 1
        if(i.ge.i_max+1) then
          write(op,100) 'No peaks were found in spectrum.'
          goto 900
        endif
* locate the first minimum after the huge peak at the origin
        if(spec(i).le.spec(i-1))
     |goto 20
      i_min = i - 1
*
* NOTE: The absolute angle is th2_low + th2_min
      th2_low = i_min * d_theta
*
  900 TRMSPC = .true.
*
      return
  100 format(1x, a)
      end
*
* ______________________________________________________________________
* Title: TST_MIR
* Author: MMJT
* Date: 30 July 1989; 22 Feb 1995
* Identifies the presence of mirror planes which contain the streaky
* axis. The integer argument 'mir_sym' can have one of two values,
* 1 or 2.
* mir_sym = 1 is the plane containing the cell sides h and l.
* mir_sym = 2 is the plane containing the cell sides k and l.
* For rotational symmetry greater than 2, one mirror implies the
* presence of the other. For diffraction symmetry group No 3 (2/M),
* only one or the other can occur, but we must test for both.
* TST_MIR returns '.true.' if a mirror was found, '.false.' otherwise.
*
*      ARGUMENTS:
*            mir_sym  -  Plane across which we wish to test for mirror
*                        symmetry. Takes the values 1 or 2. (input).
*            idum     -  parameter used by RAN3. Is -ve if RAN3
*                        is to be reset. (input).
*            ok       -  logical flag indicating all went well.
*                                                      (output).
*
*      COMMON VARIABLES:
*            uses:  a0, b0, c0, d0, lambda, DoSymDump, no_trials,
*                   max_angle, PI, PI2, RAD2DEG, cell_gamma, check_sym
*                   h_bnd, k_bnd, l_bnd, tolerance, tiny_inty
*
*        modifies:  max_var
*
*      TST_MIR returns logical .true. if the diffraction intensities
*      have the mirror symmetry about the plane requested.
* ______________________________________________________________________
*
      logical function TST_MIR(mir_sym, idum, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 mir_sym, idum
      logical ok
*
      logical is_good, match, eq_sides
      integer*4 i, h, k, h_tmp, k_tmp
      logical cell90, cell120
      real*8 RAN3, PNTINT, S, ANGLE
      real*8 tiny, l
      real*8 i_avg, tol, i1, i2, variance, rel_var
      parameter (tiny = FIVE * eps4)
*
* external functions
      external RAN3, PNTINT
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
* ANGLE is the Bragg angle (in radians) of the h,k,l plane
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
*
      cell90  = abs(cell_gamma - HALF*PI) .lt. HALF*PI*tiny
      cell120 = abs(cell_gamma - PI2/THREE) .lt. PI2*tiny/THREE
      eq_sides = abs(cell_a - cell_b) .le. HALF*eps6*(cell_a + cell_b)
      TST_MIR = .false.
      is_good = .false.
      if(mir_sym.lt.1 .or. mir_sym.gt.3) goto 900
      if(DoSymDump) then
        write(sy,200)
        if(.not.cell90 .and. .not.cell120) then
          write(sy,230) 'cell angle = ', cell_gamma * RAD2DEG,
     |      ' degrees. NO HORIZONTAL MIRRORS ARE LIKELY.'
          goto 900
        endif
      endif
*
      if(DoSymDump) then
        if(mir_sym.eq.1) then
          write(sy,199) 'Testing for mirror about the h-l plane'
        else if(mir_sym.eq.2) then
          write(sy,199) 'Testing for mirror about the k-l plane'
        else if(mir_sym.eq.3) then
          write(sy,199) 'Testing for mirror about the h=k,l plane'
        endif
        write(sy,200)
        write(sy,210)
      endif
      do 10 i = 1, no_trials
* get usable h,k >= 0
   20   if(mir_sym.eq.1) then
          h_tmp = int( dble(h_bnd + 1) * RAN3(idum) )
   30     k_tmp = int( dble(k_bnd + 1) * RAN3(idum) )
          if(k_tmp.eq.0) goto 30
        else if(mir_sym.eq.2) then
          k_tmp = int( dble(k_bnd + 1) * RAN3(idum) )
   40     h_tmp = int( dble(h_bnd + 1) * RAN3(idum) )
          if(h_tmp.eq.0) goto 40
        else if(mir_sym.eq.3) then
          k_tmp = int( dble(k_bnd + 1) * RAN3(idum) )
   45     h_tmp = int( dble(h_bnd + 1) * RAN3(idum) )
          if(h_tmp.eq.k_tmp) goto 45
        endif
* get usable l > 0
   50   l = l_bnd * RAN3(idum)
        if(abs(l).le.eps2) goto 50
* make sure we are not sampling at too high an angle
        if(TWO * ANGLE(h_tmp, k_tmp, l) .gt. max_angle) goto 20
* I(h,k,l)
        h = h_tmp
        k = k_tmp
        i1 = PNTINT(h, k, l, ok)
        if(.not.ok) goto 999
        if(DoSymDump) write(sy,220) h, k, l, i1
* mirror on h-l plane
        if(mir_sym.eq.1) then
* is the cell angle equal to 90 degrees
          if(cell90) then
* I(h,-k,l), rectangular cell
            h =  h_tmp
            k = -k_tmp
            i2 = PNTINT(h, k, l, ok)
            if(.not.ok) goto 999
            if(DoSymDump) write(sy,220) h, k, l, i2
          else if(cell120) then
* I(h+k,-k,l), hexagonal cell
            h =  h_tmp + k_tmp
            k = -k_tmp
            i2 = PNTINT(h, k, l, ok)
            if(.not.ok) goto 999
            if(DoSymDump) write(sy,220) h, k, l, i2
          endif
* else mirror on k-l plane, mir = 2
        else if(mir_sym.eq.2) then
* is the cell angle equal to 90 degrees
          if(cell90) then
* I(-h,k,l), rectangular cell
            h = -h_tmp
            k =  k_tmp
            i2 = PNTINT(h, k, l, ok)
            if(.not.ok) goto 999
            if(DoSymDump) write(sy,220) h, k, l, i2
          else if(cell120) then
* I(-h,h+k,l), hexagonal cell
            h = -h_tmp
            k =  h_tmp + k_tmp
            i2 = PNTINT(h, k, l, ok)
            if(.not.ok) goto 999
            if(DoSymDump) write(sy,220) h, k, l, i2
          endif
* else mirror on hk-l plane, mir = 3
        else if(mir_sym.eq.3) then
* The following if block is redundant, and in special
* cases fails to print to the .sym file. mmjt 3/18/04
* is the cell square
*          if(cell90 .and. eq_sides) then
* I(-h,k,l), square cell
*            h = k_tmp
*            k = h_tmp
*            i2 = PNTINT(h, k, l, ok)
*            if(.not.ok) goto 999
*            if(DoSymDump) write(sy,220) h, k, l, i2
*          else if(cell120) then
* I(-h,h+k,l), hexagonal cell
*            h = k_tmp
*            k = h_tmp
*            i2 = PNTINT(h, k, l, ok)
*            if(.not.ok) goto 999
*            if(DoSymDump) write(sy,220) h, k, l, i2
*          endif
          h = k_tmp
          k = h_tmp
          i2 = PNTINT(h, k, l, ok)
          if(.not.ok) goto 999
          if(DoSymDump) write(sy,220) h, k, l, i2
        endif
* compare mirrored intensities
        i_avg = HALF * (i1 + i2)
        variance = HALF * (abs(i_avg-i1) + abs(i_avg-i2))
* Be careful intensities are not actually zero
        if(i_avg.lt.tiny_inty) then
          tol = tiny_inty
        else
          tol = i_avg * tolerance
          rel_var = variance / i_avg
          if(rel_var.gt.max_var) max_var = rel_var
        endif
        match = abs(i_avg-i1).lt.tol .and. abs(i_avg-i2).lt.tol
        is_good = (i.eq.1 .or. is_good) .and. match
        if(DoSymDump) then
          write(sy,270)
          write(sy,240) i_avg
          write(sy,360) variance, HUNDRED * variance / i_avg
*          if(.not.check_sym) then
            write(sy,260) tol
            if(match) then
              if(mir_sym.eq.1) then
                write(sy,250)
              else if(mir_sym.eq.2) then
                write(sy,280)
              else if(mir_sym.eq.3) then
                write(sy,285)
              endif
            else
              if(mir_sym.eq.1) then
                write(sy,290)
              else if(mir_sym.eq.2) then
                write(sy,300)
              else if(mir_sym.eq.3) then
                write(sy,305)
              endif
            endif
*          endif
          write(sy,200)
        endif
   10 continue
      TST_MIR = is_good
*
      if(DoSymDump) then
*        if(.not.check_sym) then
          if(mir_sym.eq.1) then
            if(is_good) then
              write(sy,310)
            else
              write(sy,320)
            endif
          else if(mir_sym.eq.2) then
            if(is_good) then
              write(sy,330)
            else
              write(sy,340)
            endif
          else if(mir_sym.eq.3) then
            if(is_good) then
              write(sy,345)
            else
              write(sy,346)
            endif
          endif
*        endif
        write(sy,200)
      endif
*
  900 return
  999 write(op,199) 'ERROR in intensity calculation in TST_MIR'
      write(op,350) '   at h,k,l = ', h,',', k,',', l
      return
  199 format(1x, a)
  200 format(' ')
  210 format(1x, '  h', 5x, 'k', 7x, 'l', 20x, 'Intensity')
  220 format(1x, i3, 3x, i3, 2x, f9.4, 5x, f22.6)
  230 format(1x, a, f6.2, a)
  240 format(6x, 'Average Intensity = ', f22.6)
  250 format(1x, 'Intensities are consistent with an h-l mirror plane')
  260 format(1x, 'Intensity tolerance = +/-', f22.6)
  270 format(26x, '----------------------')
  280 format(1x, 'Intensities are consistent with a k-l mirror plane')
  285 format(1x, 'Intensities are consistent with a h=k,l mirror plane')
  290 format(1x, 'Intensities not consistent with an h-l mirror plane')
  300 format(1x, 'Intensities not consistent with a k-l mirror plane')
  305 format(1x, 'Intensities not consistent with a h=k,l mirror plane')
  310 format(1x, 'THERE IS A MIRROR ABOUT THE H-L PLANE')
  320 format(1x, 'THERE IS NO MIRROR ABOUT THE H-L PLANE')
  330 format(1x, 'THERE IS A MIRROR ABOUT THE K-L PLANE')
  340 format(1x, 'THERE IS NO MIRROR ABOUT THE K-L PLANE')
  345 format(1x, 'THERE IS A MIRROR ABOUT THE H=K,L PLANE')
  346 format(1x, 'THERE IS NO MIRROR ABOUT THE H=K,L PLANE')
  350 format(1x, a, i3, a, i3, a, f7.2)
  360 format(1x,'  Average variation = +/-', f22.6,'  (+/-',g9.2,'%)')
      end
*
* ______________________________________________________________________
* Title: TST_ROT
* Author: MMJT
* Date: 30 July 1989; 22 Feb 1995
* Identifies the rotational symmetry of the diffraction pattern.
* The rotational symmetry to test for is passed as 'rot_sym'.
* The values accepted for 'rot_sym' are 2, 3 and 4. If the diffraction
* intensity has the symmetry requested, TST_ROT returns '.true.'. If
* the pattern does not have the requested symmetry, or if an illegal
* value was passed (i.e. rot_sym = 5), TST_ROT returns '.false.'.
* NOTE. A 6-fold axis is found by calling TST_ROT twice: once for
* rot_sym = 2 and again for rot_sym = 3.
*
*      ARGUMENTS:
*            rot_sym  -  Rotational symmetry to test for.
*                        Accepts the values 2, 3 or 4. (input).
*            idum     -  parameter used by RAN3. Is -ve if RAN3
*                        is to be reset. (input).
*            ok       -  logical flag indicating all went well.
*                                                      (output).
*
*      COMMON VARIABLES:
*            uses:  a0, b0, c0, d0, lambda, DoSymDump, no_trials,
*                   max_angle, check_sym, h_bnd, k_bnd, l_bnd
*                   tolerance
*
*        modifies:  max_var
*
*      TST_ROT returns logical .true. if the diffraction intensities
*      have the requested rotational symmetry.
* ______________________________________________________________________
*
      logical function TST_ROT(rot_sym, idum, ok)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 rot_sym, idum
      logical ok
*
      logical is_good, match
      integer*4 i, h, k, h_tmp, k_tmp
      real*8 RAN3, PNTINT, S, ANGLE
      real*8 l, i_avg, tol
      real*8 i1, i2, i3, i4, variance, rel_var
*
* external functions
      external RAN3, PNTINT
*
* statement functions
* S is the value of 1/d**2 at hkl
      S(h,k,l) = h*h*a0 + k*k*b0 + l*l*c0 + h*k*d0
      ANGLE(h,k,l) = asin(HALF * lambda * sqrt(S(h,k,l)))
*
      TST_ROT = .false.
      is_good = .false.
* Is rot valid?
      if(rot_sym.lt.2 .or. rot_sym.gt.4) goto 900
* Now test for rotational symmetry.
* 2-fold and 4-fold
* Avoid both h, k = 0. Also, avoid l = 0, since Friedel's law will
* create a pseudo 2-fold.
      if(rot_sym.eq.2 .or. rot_sym.eq.4) then
        if(DoSymDump) then
          write(sy,210)
          write(sy,330) 'Testing for ', rot_sym, '-fold axis'
          write(sy,210)
        endif
        do 10 i = 1, no_trials
          if(DoSymDump) write(sy,220)
* get usable h,k >= 0
   20     h_tmp = int( dble(h_bnd + 1) * RAN3(idum) )
          k_tmp = int( dble(k_bnd + 1) * RAN3(idum) )
          if(h_tmp.eq.0 .and. k_tmp.eq.0) goto 20
* get usable l > 0
   30     l = l_bnd * RAN3(idum)
* keep l off the l = 0 plane, else we might confuse the inversion
* with a 2-fold
          if(abs(l).le.eps2) goto 30
* make sure we are not sampling at too high an angle
          if(TWO * ANGLE(h_tmp, k_tmp, l) .gt. max_angle) goto 20
* I(h,k,l)
          h = h_tmp
          k = k_tmp
          i1 = PNTINT(h, k, l, ok)
          if(.not.ok) goto 999
          if(DoSymDump) write(sy,230) h, k, l, i1
* I(-h,-k,l)
          h = -h_tmp
          k = -k_tmp
          i2 = PNTINT(h, k, l, ok)
          if(.not.ok) goto 999
          if(DoSymDump) write(sy,230) h, k, l, i2
* compare 2-fold intensities
          if(rot_sym.eq.2) then
            i_avg = HALF * (i1 + i2)
            variance = HALF * (abs(i_avg-i1) + abs(i_avg-i2))
* Be careful intensities are not actually zero
            if(i_avg.lt.tiny_inty) then
              tol = tiny_inty
            else
              tol = i_avg * tolerance
              rel_var = variance / i_avg
              if(rel_var.gt.max_var) max_var = rel_var
            endif
            match = abs(i_avg-i1).lt.tol .and. abs(i_avg-i2).lt.tol
            is_good = (i.eq.1 .or. is_good) .and. match
          else
* I(-k,h,l)
            h = -k_tmp
            k =  h_tmp
            i3 = PNTINT(h, k, l, ok)
            if(.not.ok) goto 999
            if(DoSymDump) write(sy,230) h, k, l, i3
* I(k,-h,l)
            h =  k_tmp
            k = -h_tmp
            i4 = PNTINT(h, k, l, ok)
            if(.not.ok) goto 999
            if(DoSymDump) write(sy,230) h, k, l, i4
* compare 4-fold intensities
            i_avg = QUARTER * (i1 + i2 + i3 + i4)
            variance = QUARTER * (abs(i_avg-i1) + abs(i_avg-i2) +
     |                           abs(i_avg-i3) + abs(i_avg-i4))
* Be careful intensities are not actually zero
            if(i_avg.lt.tiny_inty) then
              tol = tiny_inty
            else
              tol = i_avg * tolerance
              rel_var = variance / i_avg
              if(rel_var.gt.max_var) max_var = rel_var
            endif
            match = abs(i_avg-i1).lt.tol .and. abs(i_avg-i2).lt.tol
     |        .and. abs(i_avg-i3).lt.tol .and. abs(i_avg-i4).lt.tol
            is_good = (i.eq.1.or.is_good) .and. match
          endif
*
          if(DoSymDump) then
            write(sy,240)
            write(sy,250) i_avg
            write(sy,260) variance, HUNDRED * variance / i_avg
*            if(.not.check_sym) then
              write(sy,270) tol
              if(match) then
                write(sy,280) rot_sym
              else
                write(sy,290) rot_sym
              endif
*            endif
            write(sy,210)
          endif
   10   continue
        TST_ROT = is_good
        goto 900
      endif
* 3-fold
* Avoid both h, k = 0.
      if(rot_sym.eq.3) then
        if(DoSymDump) then
          write(sy,200) rot_sym
          write(sy,210)
          write(sy,220)
        endif
        do 40 i = 1, no_trials
* get usable h,k >= 0
   50     h_tmp = int( dble(h_bnd + 1) * RAN3(idum) )
          k_tmp = int( dble(k_bnd + 1) * RAN3(idum) )
          if(h_tmp.eq.0 .and. k_tmp.eq.0) goto 50
* get l (l=0 is allowed)
          l = l_bnd * RAN3(idum)
* make sure we are not sampling at too high an angle
          if(TWO * ANGLE(h_tmp, k_tmp, l) .gt. max_angle) goto 50
* I(h,k,l)
          h = h_tmp
          k = k_tmp
          i1 = PNTINT(h, k, l, ok)
          if(.not.ok) goto 999
          if(DoSymDump) write(sy,230) h, k, l, i1
* I(-h-k,h,l)
          h = -(h_tmp + k_tmp)
          k = h_tmp
          i2 = PNTINT(h, k, l, ok)
          if(.not.ok) goto 999
          if(DoSymDump) write(sy,230) h, k, l, i2
* I(k,-h-k,l)
          h = k_tmp
          k = -(h_tmp + k_tmp)
          i3 = PNTINT(h, k, l, ok)
          if(DoSymDump) write(sy,230) h, k, l, i3
* compare intensities
          i_avg = (i1 + i2 + i3) / THREE
          variance =
     |      (abs(i_avg-i1) + abs(i_avg-i2) + abs(i_avg-i3)) / THREE
* Be careful intensities are not actually zero
          if(i_avg.lt.tiny_inty) then
            tol = tiny_inty
          else
            tol = i_avg * tolerance
            rel_var = variance / i_avg
            if(rel_var.gt.max_var) max_var = rel_var
          endif
          match = abs(i_avg-i1).lt.tol .and. abs(i_avg-i2).lt.tol
     |      .and. abs(i_avg-i3).lt.tol
          is_good = (i.eq.1 .or. is_good) .and. match
          if(DoSymDump) then
            write(sy,240)
            write(sy,250) i_avg
            write(sy,260) variance, HUNDRED * variance / i_avg
*            if(.not.check_sym) then
              write(sy,270) tol
              if(match) then
                write(sy,280) rot_sym
              else
                write(sy,290) rot_sym
              endif
*            endif
            write(sy,210)
          endif
   40   continue
        TST_ROT = is_good
      endif
*
  900 if(DoSymDump) then
*        if(.not.check_sym) then
          if(is_good) then
            write(sy,300) rot_sym
          else
            write(sy,310) rot_sym
          endif
*        endif
        write(sy,210)
      endif
      return
*
  999 write(op,400) 'ERROR in intensity calculation in TST_ROT'
      write(op,320) '   at h,k,l = ', h,',', k,',', l
      return
  200 format(1x, 'Testing for a ', i1, '-fold axis')
  210 format(' ')
  220 format(1x, '  h', 5x, 'k', 7x, 'l', 20x, 'Intensity')
  230 format(1x, i3, 3x, i3, 2x, f9.4, 5x, f22.6)
  240 format(26x, '----------------------')
  250 format(6x, 'Average Intensity = ', f22.6)
  260 format(1x, '  Average variation = +/-',f22.6,'  (+/-',g9.2,'%)')
  270 format(1x, 'Intensity tolerance = +/-', f22.6)
  280 format(1x, 'Intensities are consistent with a ', i1, '-fold')
  290 format(1x, 'Intensities are not consistent with a ', i1, '-fold')
  300 format(1x, 'INTENSITY DISTRIBUTION HAS A ', i1, '-FOLD AXIS')
  310 format(1x, 'INTENSITY DISTRIBUTION HAS NO ', i1, '-FOLD AXIS')
  320 format(1x, a, i3, a, i3, a, f7.2)
  330 format(1x, a, i1, a)
  400 format(1x, a)
      end
*
* ______________________________________________________________________
** Title: WRTLNE
** Author: MMJT
** Date: 4 Oct 1989
** Description: This subroutine writes out the contents of
** the character string 'line' with an appropriate preamble.
** This routine strips off trailing blanks.
**
**      ARGUMENTS:
**            line     -  Line of characters to write out. (input).
** ______________________________________________________________________
**
*      subroutine WRTLNE(line)
*      include 'DIFFaX.par'
**     save
**
*      character*(*) line
**
*      integer*4 lin_len, i
**
*      lin_len = len(line)
*      i = lin_len
*   10 if(i.gt.0) then
*        if(line(i:i).eq.' ') then
*          if(i.gt.1) then
*            i = i - 1
*            goto 10
*          endif
*        endif
*      endif
**
*      write(op,100) 'The data on the current line was read as:'
*      write(op,101) '  ''', line(1:i), ''''
**
*      return
*  100 format(1x, a)
*  101 format(1x, 3a)
*      end
**
** ______________________________________________________________________
** Title: WRTSAD
** Author: MMJT
** Date: 29 Oct 1989; 1st Mar 2000; 19th May 2010
** Description: This subroutine writes the selected area diffraction
** pattern (SADP) data to a binary file.
** Binary file output has been problematical because some compilers add
** record-length specifiers at the beginning and ending of each line.
** Sometimes the specifiers are 2-byte, or 4-byte, and even 0-byte!
** A statement to screen specifies how to compute the true row width.
**
**      ARGUMENTS:
**            outfile  -  The name of the file that the binary
**                        selected area diffraction data is to be
**                        written to. (input).
**            view     -  Choice of beam direction. (input).
**                              1  =  normal to the plane k = 0
**                              2  =  normal to the plane h = 0
**                              3  =  normal to the plane h = k
**                              4  =  normal to the plane h = -k
**            l_upper  -  Upper limit of l. (input).
**            hk_lim   -  Upper limit of h (or k). (input).
**            ok       -  logical flag indicating all went well.
**                                                      (output).
**
**      COMMON VARIABLES:
**            uses:  a0, b0, c0, d0, sadblock, scaleint, has_l_mirror
**
**        modifies:  spec
** ______________________________________________________________________
**
*      subroutine WRTSAD(outfile, view, l_upper, hk_lim, ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      character*(*) outfile
*      integer*4 view, hk_lim
*      real*8 l_upper
*      logical ok
**
*      integer*4 i, j, p1, p2, LENGTH
*      real*8 rowint(SADSIZE), incr, sigma, x
**
** external function
*      external LENGTH
**
** external subroutine (Some compilers need them declared external)
**      external SMUDGE
**
*      write(op,102) 'Writing SADP data to file ''',
*     |  outfile(1:LENGTH(outfile)), '''. . .'
**
** set standard deviation
*      sigma = ZERO
**
** Establish scaling in reciprocal space, the number of pixels per unit
** first get number of pixels per l
*      incr = dble(SADSIZE/2) / l_upper
*      if(view.eq.1) then
** number of pixels per unit h
*        incr = incr * sqrt(a0 / c0)
*      else if(view.eq.2) then
** number of pixels per unit k
*        incr = incr * sqrt(b0 / c0)
*      else if(view.eq.3) then
** number of pixels per unit along (h = k)
*        incr = incr * sqrt((a0 + b0 + d0) / c0)
*      else if(view.eq.4) then
** number of pixels per unit along (h = -k)
*        incr = incr * sqrt((a0 + b0 - d0) / c0)
*      endif
**
*      do 10 j = sadblock - 1, 0, -1
*        do 20 i = 1, SADSIZE
*          rowint(i) = ZERO
*   20   continue
*        do 30 i = 0, hk_lim
*          p1 = SADSIZE/2 + nint(i*incr) + 1
*          p2 = SADSIZE/2 - nint(i*incr) + 1
**
** cycle if we have gone out of bounds
*          if(p1.gt.SADSIZE .or. p1.lt.0) then
*            goto 30
*          endif
*          if(p2.gt.SADSIZE .or. p2.lt.0) then
*            goto 30
*          endif
**
*          x = spec(i*sadblock + j + 1) * scaleint
** handle saturated pixels
*          if(x.gt.maxsad) x = maxsad
*          rowint(p1) = x
*          if(has_l_mirror) then
*            rowint(p2) = x
*          else
*            x = spec(i*sadblock + sadblock - j) * scaleint
** handle saturated pixels
*            if(x.gt.maxsad) x = maxsad
*            rowint(p2) = x
*          endif
*   30   continue
*        call SMUDGE(rowint, SADSIZE, sigma, ok)
*        if(.not.ok) goto 999
*        if(bitdepth.eq.8) then
*          write(sa,err=900) (char(nint(rowint(i))), i = 1, SADSIZE)
*        else
** Known bug: Writing 16-bit binary is problematic on Intel processors
*            write(sa,err=900) 
*     | (char(int(rowint(i)/256)),
*     |  char(int(rowint(i)-(int(rowint(i)/256))*256)),i=1,SADSIZE)
**            write(sa,err=900) (iidnnt(rowint(i)), i = 1, SADSIZE)
*        endif
*   10 continue
**
** Repeat, in reverse for the bottom half if data had a mirror.
*      if(has_l_mirror) then
*        do 40 j = 1, sadblock - 1
*          do 50 i = 1, SADSIZE
*            rowint(i) = ZERO
*   50     continue
*          do 60 i = 0, hk_lim
*            p1 = SADSIZE/2 + nint(i*incr) + 1
*            p2 = SADSIZE/2 - nint(i*incr) + 1
**
** cycle if we have gone out of bounds
*            if(p1.gt.SADSIZE .or. p1.lt.0) then
*              goto 60
*            endif
*            if(p2.gt.SADSIZE .or. p2.lt.0) then
*              goto 60
*            endif
**
*            x = spec(i*sadblock + j + 1) * scaleint
** handle saturated pixels
*            if(x.gt.maxsad) x = maxsad
*            rowint(p1) = x
*            rowint(p2) = x
*   60     continue
*          call SMUDGE(rowint, SADSIZE, sigma, ok)
*          if(.not.ok) goto 999
*          if(bitdepth.eq.8) then
*            write(sa,err=900) (char(nint(rowint(i))), i = 1, SADSIZE)
*          else
** Known bug: Writing 16-bit binary is problematic on Intel processors
*            write(sa,err=900) 
*     | (char(int(rowint(i)/256)),
*     |  char(int(rowint(i)-(int(rowint(i)/256))*256)),i=1,SADSIZE)
**            write(sa,err=900) (iidnnt(rowint(i)), i = 1, SADSIZE)
*          endif
*   40   continue
** write a blank last line to make the array SADSIZE x SADSIZE square
*        if(bitdepth.eq.8) then
*          write(sa,err=900) (char(0), i = 1, SADSIZE)
*        else
*          write(sa,err=900) (char(0), i = 1, 2*SADSIZE)
*        endif
*      endif
**
*      if(sa.ne.op) close(unit = sa, err = 990)
**
** MMJT 19th May 2010
*      write(op,103) SADSIZE, ' x ', SADSIZE, ' pixels: ',
*     |              bitdepth, ' bits deep.'
*      write(op,100) ' (Caution: Because of a compiler-dependent issue'
*      write(op,100) '  the true image width in pixels is found by'
*      if(bitdepth.eq.8) then
*        write(op,100) '  dividing the file size in bytes by 256.)'
*      else
*        write(op,100) '  dividing the file size in bytes by 512.'
*        write(op,100) '  The byte-order is big-endian.)'
*      endif
**
*      return
*  900 write(op,100) 'ERROR: problems writing to binary SADP file.'
*      ok = .false.
*      return
*  990 write(op,100) 'ERROR: problems closing binary SADP file.'
*      ok = .false.
*      return
*  999 write(op,100) 'ERROR: SMUDGE returned error to WRTSAD.'
*      write(op,101) i, j
*      return
*  100 format(1x, a)
*  101 format(5x, 'with local variables i = ', i5, ' j = ', i5)
*  102 format(1x, 3a)
*  103 format(1x, i4, a, i4, a, i2, a)
*      end
** ______________________________________________________________________
** Title: WRTSPC
** Author: MWD and MMJT
** Date: 18 Aug 1988; 7 Mar 1995
** Description:  This routine writes the spectrum arrays to file.
**
**      ARGUMENTS:
**            spcfile  -  The name of the output data file. (input).
**            ok       -  logical flag indicating all went well.
**                                                      (output).
**
**      COMMON VARIABLES:
**            uses:  th2_min, th2_max, d_theta, spec, brd_spc, blurring
**                   NONE, RAD2DEG
**
**        modifies:  no COMMON variables are modified
** ______________________________________________________________________
**
*      subroutine WRTSPC(spcfile, ok)
*      include 'DIFFaX.par'
*      include 'DIFFaX.inc'
**
*      character*(*) spcfile
*      logical ok
**
*      character*1 tab
*      integer*4 i, n_low, n_high, LENGTH
*      real*8 theta
**
** external function
*      external LENGTH
**
*      n_low  = 1
*      n_high = int(HALF*(th2_max-th2_min)/d_theta) + 1
**
*      tab = char(9)
*      theta = th2_min * RAD2DEG
*      write(op,200) 'Writing spectrum data to file ''',
*     |      spcfile(1:LENGTH(spcfile)), '''. . .'
**      do 10 i = int(HALF*th2_min / d_theta) + 1,
**     |                  int(HALF*th2_max / d_theta) + 1
*      if(blurring.eq.NONE) then
*        do 10 i = n_low, n_high
*          write(sp,101,err=50) theta, tab, spec(i)
*          theta = theta + TWO * RAD2DEG * d_theta
*   10   continue
*      else
*        do 20 i = n_low, n_high
*          write(sp,100,err=50) theta, tab, spec(i), tab, brd_spc(i)
*          theta = theta + TWO * RAD2DEG * d_theta
*   20   continue
*      endif
**
*      if(sp.ne.op) close(sp,err=60)
*      write(op,202) 'Spectrum written.'
*      write(op,202)
*      return
*   50 write(op,202) 'Unable to write to spectrum file.'
*      if(sp.ne.op) close(sp,err=60)
*      ok = .false.
*      return
*   60 write(op,202) 'Unable to close spectrum file.'
*      ok = .false.
*      return
*  100 format(1x, e12.5, 2(a, g13.6))
*  101 format(1x, e12.5, a, g13.6)
*  200 format(1x, 3a)
*  201 format(1x, a, i2, 2a)
*  202 format(1x, a)
*      end
**
** ______________________________________________________________________
* Title: XYPHSE
* Author: MMJT
* Date: 8 June 1990
* Description:  This routine pre-calculates the h and k components of
* phases for each atom. This is called when l = 0, since h and k
* are held constant while l is varied along each streak.
*
*      ARGUMENTS:
*            h  -  reciprocal lattice vector h-component. (input).
*            k  -  reciprocal lattice vector k-component. (input).
*
*      COMMON VARIABLES:
*            uses:  n_actual, l_n_atoms, a_pos, l_actual, n_layers
*
*        modifies:  hx_ky
* ______________________________________________________________________
*
      subroutine XYPHSE(h, k)
      include 'DIFFaX.par'
      include 'DIFFaX.inc'
*
      integer*4 h, k
*
      integer*4 i, m
*
      do 10 m = 1, n_actual
        do 20 i = 1, l_n_atoms(m)
          hx_ky(i,m) = h*a_pos(1,i,m) + k*a_pos(2,i,m)
   20   continue
   10 continue
*
      return
      end
*
* ______________________________________________________________________
* Title: YRDSTK
* Author: MMJT
* Date: 12 Aug 1989
* Description: YRDSTK checks that y is a multiple of x.
*
*      ARGUMENTS:
*            x   -  reference value. (input).
*            y   -  value to be tested. (input)
*            ok  -  logical flag set to .false. if y is zero. (output).
*
*      YRDSTK returns logical .true. if y is a multiple of x.
* ______________________________________________________________________
*
      logical function YRDSTK(x, y, ok)
      include 'DIFFaX.par'
*     save
*
      real*8 x, y
      logical ok
*
      real*8 tmp
*
      YRDSTK = .false.
      if(y.eq.ZERO) goto 999
      if(x.eq.ZERO) then
* This is the first visit to YRDSTK. Set things up.
        x = y
        YRDSTK = .true.
      else
        tmp = y / x
        if(abs(nint(tmp)-tmp) .le. eps3*tmp) YRDSTK = .true.
      endif
*
      return
  999 ok = .false.
      return
      end
*
