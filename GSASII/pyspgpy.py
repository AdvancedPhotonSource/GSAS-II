"""
pyspgpy: all-Python space group interpreter
===========================================

A Python re-implementation of A.C. Larson's Fortran space-group routines
(pyspg.for, sgroupnp.for, sglatc.for, sglcen.for, sglpak.for, sgmtml.for,
sgoprn.for, sgrmat.for and sgtrcf.for, National Research Council of
Canada), which generate the symmetry operators of a space group from its
Hermann-Mauguin symbol.

Translated from Fortran by Copilot/Kimi-K3 in a few steps (BHT, 11/26)

The entry point is :func:`sgforpy`, which takes a space-group symbol (a
character string) and returns the list
``[LAUE, SGINV, SGLATT, SGUNIQ, SGPOL, SGNOPS, SGMTRX, SGTRNS, SGGEN, IERR]``:

* LAUE:   Laue group number (1=-1, 2=2/m, 3=mmm, 4=4/m, 5=4/mmm, 6=R3R,
          7=R3MR, 8=3, 9=3m1, 10=31m, 11=6/m, 12=6/mmm, 13=m3, 14=m3m)
* SGINV:  centrosymmetry flag (0 = acentric, 1 = centric)
* SGLATT: lattice centering (1=P, 2=A, 3=B, 4=C, 5=I, 6=F, 7=R)
* SGUNIQ: unique axis for monoclinic groups (1=a, 2=b, 3=c; -1 for
          rhombohedral in the hexagonal setting; 4 on error exit)
* SGPOL:  polar axis flag (1=x, 2=y, 3=xy, 4=z, 5=xz, 6=yz, 7=xyz, 8=111)
* SGNOPS: number of symmetry operators generated
* SGMTRX: (SGNOPS,3,3) array of rotation matrices
* SGTRNS: (SGNOPS,3) array of translation vectors
* SGGEN:  (SGNOPS,) array of symmetry operation flags
* IERR:   error flag (0 = no error)

Internally the operators are represented exactly: rotation matrices as
integers and translations as integer twelfths, so no floating-point
rounding occurs.
"""
import numpy as np


class _SgError(Exception):
    """Bail out to the error exit carrying the Fortran IER error code."""


# --- symbol parsing ----------------------------------------------------------

# characters that may appear in a symbol; tab is treated as a space
_LEGAL = set('CBAPFIRMND123456-/H. ')
# the glide/mirror letters; a letter terminates its field
_LETTERS = 'CBAPFIRMND'
# lower case letters are folded to upper case
_FOLD = {c: c.upper() for c in 'cbapfirmndh'}
# lattice character -> centering flag (1=P, 2=A, 3=B, 4=C, 5=I, 6=F, 7=R)
_LATTICE = {'P': 1, 'A': 2, 'B': 3, 'C': 4, 'I': 5, 'F': 6, 'R': 7}
# the numeric character codes of the original Fortran, for a few relational
# tests ('' is the padding token)
_CODE = {'': 0}
_CODE.update({c: i for i, c in enumerate(' CBAPFIRMND123456-/H.', 1)})
# operator character -> operator kind (the Fortran I = ABS(code-5) dispatch):
# 1-5: A,B,C,M,N mirrors/axial glides (quirk: F,I,R act as A,B,C);
# 6: d-glide; 7-12: 1-,2-,3-,4-,5-,6-fold axes; 13-15: '-', '/', 'H'
# (which generate no operator)
_OP = {'A': 1, 'F': 1, 'B': 2, 'I': 2, 'C': 3, 'R': 3, 'M': 4, 'N': 5,
       'D': 6, '1': 7, '2': 8, '3': 9, '4': 10, '5': 11, '6': 12,
       '-': 13, '/': 14, 'H': 15}


def _parse_symbol(spg):
    """Break the space group symbol into fields of operator tokens.

    Returns (fields, k): fields is a list of four fields, each a list of
    four tokens with '' padding, and k is the number of fields found.
    """
    fields = []
    field = []
    for ch in spg[:20].replace('\t', ' '):
        ch = _FOLD.get(ch, ch)
        if ch == ' ':
            if field:                           # squeeze out extra spaces
                fields.append(field)
                field = []
        elif ch in _LEGAL:
            field.append(ch)
            if ch in _LETTERS or len(field) == 4:
                fields.append(field)            # a letter or a 4th token ends a field
                field = []
        else:
            raise _SgError(29)
        if len(fields) == 4:                    # at most 4 fields are used
            break
    if field and len(fields) < 4:
        fields.append(field)
    k = len(fields)
    fields += [[]] * (4 - k)
    return [f + [''] * (4 - len(f)) for f in fields], k


# --- symmetry operators ------------------------------------------------------

# lattice centering vectors in twelfths (1=P, 2=A, 3=B, 4=C, 5=I, 6=F, 7=R)
_CENTERING = {
    1: [(0, 0, 0)],
    2: [(0, 0, 0), (0, 6, 6)],
    3: [(0, 0, 0), (6, 0, 6)],
    4: [(0, 0, 0), (6, 6, 0)],
    5: [(0, 0, 0), (6, 6, 6)],
    6: [(0, 0, 0), (0, 6, 6), (6, 0, 6), (6, 6, 0)],
    7: [(0, 0, 0), (4, 8, 8), (8, 4, 4)],
}

# rotation matrices used to seed the group
_IDENTITY = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
_ZERO = ((0, 0, 0), (0, 0, 0), (0, 0, 0))
_ZXY = ((0, 0, 1), (1, 0, 0), (0, 1, 0))        # 3-fold along [111]: z,x,y
_YZX = ((0, 1, 0), (0, 0, 1), (1, 0, 0))        # 3-fold along [111]: y,z,x
_M_A = ((-1, 0, 0), (0, 1, 0), (0, 0, 1))       # mirror normal to A
_M_B = ((1, 0, 0), (0, -1, 0), (0, 0, 1))       # mirror normal to B
_M_C = ((1, 0, 0), (0, 1, 0), (0, 0, -1))       # mirror normal to C
_M_HEX = ((-1, 1, 0), (0, 1, 0), (0, 0, 1))     # mirror normal to [100], hex cell
_M_DIAG = ((0, 1, 0), (1, 0, 0), (0, 0, 1))     # diagonal mirror normal to [-110]
_2_A = ((1, 0, 0), (0, -1, 0), (0, 0, -1))      # 2-fold along A
_2_B = ((-1, 0, 0), (0, 1, 0), (0, 0, -1))      # 2-fold along B
_2_C = ((-1, 0, 0), (0, -1, 0), (0, 0, 1))      # 2-fold along C
_2_110 = ((0, 1, 0), (1, 0, 0), (0, 0, -1))     # 2-axis along [110]
_2_1M10 = ((0, -1, 0), (-1, 0, 0), (0, 0, -1))  # 2-axis along [1-10]
_2_100H = ((1, -1, 0), (0, -1, 0), (0, 0, -1))  # 2-axis along [100], hex cell
_2_210 = ((1, 0, 0), (1, -1, 0), (0, 0, -1))    # 2-axis along [210]
_3FOLD = ((0, -1, 0), (1, -1, 0), (0, 0, 1))    # 3-fold along C, hex cell
_4FOLD = ((0, -1, 0), (1, 0, 0), (0, 0, 1))     # 4-fold along C
_4BAR = ((0, 1, 0), (-1, 0, 0), (0, 0, -1))     # 4-bar along C
_6FOLD = ((1, -1, 0), (1, 0, 0), (0, 0, 1))     # 6-fold along C
_6BAR = ((-1, 1, 0), (-1, 0, 0), (0, 0, -1))    # 6-bar along C


def _pack(m):
    """A unique integer key for a 3x3 rotation matrix (the Fortran RT(5,1)
    packed matrix flag)."""
    return (81 * (2 * m[0, 0] + 3 * m[0, 1] + 4 * m[0, 2])
            + 9 * (2 * m[1, 0] + 3 * m[1, 1] + 4 * m[1, 2])
            + 2 * m[2, 0] + 3 * m[2, 1] + 4 * m[2, 2])


class _Op:
    """A symmetry operator x' = M x + T: M is an integer 3x3 rotation
    matrix and T the translation in twelfths (T/12 is the fractional
    translation)."""

    __slots__ = ('m', 't', 'key')

    def __init__(self, m, t=(0, 0, 0)):
        self.m = np.array(m, dtype=int).reshape(3, 3)
        self.t = np.array(t, dtype=int) % 12
        self.key = _pack(self.m)


def _compose(a, b):
    """The product operator of a applied after b (Fortran SGMTML)."""
    return _Op(a.m @ b.m, a.m @ b.t + a.t)


_IDENTITY_OP = _Op(_IDENTITY)


# spin-flip operation flags for each rotation matrix key (Fortran SGOPRN)
_OPFLAG = {
    193: 0, -193: 7, 185: 1, -185: 6, 139: 2, -139: 5, -131: 4, 131: 3,
    -257: 98, 257: 99, 265: 8, -265: 12, -299: 8, 299: 12, 353: 8, -353: 12,
    123: 8, -123: 12, 201: 8, -201: 12, -221: 16, 221: 23, 229: 16, -229: 23,
    -295: 16, 295: 23, 349: 16, -349: 23, 129: 16, -129: 23, 195: 16, -195: 23,
    345: 0, -345: 7, 281: 0, -281: 7, -309: 0, 309: 7, 205: 0, -205: 7,
    303: 0, -303: 7, -277: 0, 277: 7, -339: 0, 339: 7, -209: 0, 209: 7,
    -248: 0, 248: 7, 67: 0, -67: 7, -59: 16, 59: 23, 256: 16, -256: 23,
    112: 4, -112: 3, 157: 4, -157: 3, -149: 8, 149: 12, -104: 8, 104: 12,
}


def _translation_conflict(mf, opn, opm2, lcent, laueno):
    """Search for a translation conflict between two operators that share a
    rotation matrix, up to lattice centering (Fortran SGTRCF).  Returns the
    IER error code."""
    ncv = (1, 2, 3, 4, 5, 4, 1)[lcent - 1]
    jcv = (1, 1, 2, 3, 4, 1, 1)[lcent - 1]
    icenv = ((0, 0, 0), (0, 6, 6), (6, 0, 6), (6, 6, 0), (6, 6, 6))
    dr = opn.t - opm2.t
    # if the pair generates 1bar, compare against the identity instead
    m2z = _IDENTITY_OP if opn.key == -opm2.key else opm2
    icv = 1 - jcv
    tottr = 1
    while tottr != 0 and icv < ncv:
        icv += jcv
        irx1 = (dr[0] + icenv[icv - 1][0]) % 12
        iry1 = (dr[1] + icenv[icv - 1][1]) % 12
        irz1 = (dr[2] + icenv[icv - 1][2]) % 12
        if opn.m[2, 2] + m2z.m[2, 2] <= 0:
            irz1 = 0                        # Z is constrained in the unit cell
        if laueno <= 3 or mf != 3:
            # this operator does not operate along the face diagonal
            if opn.m[0, 0] + m2z.m[0, 0] <= 0:
                irx1 = 0
            if opn.m[1, 1] + m2z.m[1, 1] <= 0:
                iry1 = 0
        else:
            irx1 = (irx1 + iry1) % 12
            iry1 = 0
        tottr = 144 * irx1 + 12 * iry1 + irz1
    return 18 if tottr != 0 else 0


def _try_add(ops, cand, mf, lcent, laueno):
    """Check a candidate operator against the ones already generated and
    append it if it is new.  A candidate whose rotation matrix is the
    negative of an existing one defines a center of symmetry; a matching
    matrix with different translations signals a conflict.
    Returns (added, centric, ier)."""
    for other in ops:
        if cand.key == other.key:
            ier = 0
            if not np.array_equal(cand.t, other.t):
                ier = _translation_conflict(mf, cand, other, lcent, laueno)
            return False, False, ier
        if cand.key == -other.key:
            ier = 0
            if not np.array_equal(cand.t, other.t):
                ier = _translation_conflict(mf, cand, other, lcent, laueno)
            return False, True, ier
    ops.append(cand)
    return True, False, 0


# --- Laue group and origin definition (Fortran SGLATC) ------------------------

def _laue_group(fields, k, lcent, d):
    """Determine the Laue group and the origin-definition flags d (a 3x3
    list of integer twelfths, modified in place).
    Returns (laueno, naxis, ier, i209, id): the Laue group number, the
    unique axis flag, an error code, the body-diagonal 3-axis flag and the
    number of d-glides."""
    f2, f3, f4 = fields[1], fields[2], fields[3]
    laueno = naxis = ier = i209 = id_ = 0

    def tetragonal_tail(lbl):
        """Finish the 4/mmm origin flags (Fortran labels 1551..1559).
        Returns an IER error code."""
        while True:
            if lbl == 1551:
                d[0][2] = d[1][2] = 6
                return 0
            elif lbl == 1552:
                d[0][2] = 6
                return 0
            elif lbl == 1553:
                if lcent < 5:
                    lbl = 1555
                    continue
                if f4[0] == 'C':        # C along (110): place the D at Y=1/4
                    d[1][0] += 6
                if f2[3] != 'A' and f2[3] != 'D':   # not I 41/A or F 41/D
                    lbl = 1554
                    continue
                d[0][2] = 3
                if lcent == 5:
                    d[1][2] = 9
                return 0
            elif lbl == 1554:
                if f2[1] != '1':            # no 41 present
                    return 0
                if lcent == 6:              # F-centered
                    lbl = 1558
                    continue
                d[1][2] = 6                 # I 41 2 2
                lbl = 1557
                continue
            elif lbl == 1555:
                if lcent != 4:              # not C-centered
                    return 23
                if f2[2] == 'A' or f2[3] == 'A':    # an A normal to C
                    lbl = 1559
                    continue
                if d[0][0] == 0:
                    d[0][0] = 2 * d[1][1]
                if f4[0] == '2' and f4[1] == '1':   # a 21 on the diagonal
                    lbl = 1552
                    continue
                if f2[1] == '' or _CODE[f2[1]] > 15:
                    return 0
                if f4[0] != 'N':            # no N-glide normal to (110)
                    return 0
                d[0][0] -= 2 * d[1][1]
                lbl = 1552
                continue
            elif lbl == 1556:
                # account for translations due to diagonal symmetry operation
                if f3[0] == 'D' and lcent == 6:
                    d[2][0] = 3             # F 4bar D 2: want the 2 at Z=1/8
                if f4[0] == '2' and f4[1] == '1':
                    d[0][0] = 6             # * 4bar * 21: want the mirror at X=1/4
                if f4[0] == 'C' or f4[0] == 'N':
                    d[2][1] = 6             # C- or N-glide (110): 2-axis at Z=1/4
                if f4[0] == 'B' or f4[0] == 'N':
                    d[0][1] = 6             # B- or N-glide (110): set 2 at X=1/4
                if f4[0] != 'D':
                    return 0
                lbl = 1557
                continue
            elif lbl == 1557:
                if lcent == 5:
                    d[0][1] = 6
                d[2][1] = 9
                return 0
            elif lbl == 1558:
                d[0][2] = 3                 # F 41 * *
                d[1][2] = 9
                return 0
            elif lbl == 1559:
                d[0][2] = 3                 # C 4*/A * *
                d[1][2] = 3
                if f4[0] == 'B' or f4[0] == 'N':
                    d[0][0] = 6
                return 0

    def monoclinic_origin(im, ia, ic, na, nb, nc, un):
        """Origin definition for monoclinic groups (Fortran label 1430):
        im is the 0-based field index, ia/ic the glide letters and na/nb/nc
        and un 0-based axis indices (un is the unique axis)."""
        f = fields[im]
        if f[1] == '1':
            d[nb][un] = 6
        if f[2] == ia or f[2] == 'N':
            d[na][un] = 6
        if f[2] == ic or f[2] == 'N':
            d[nc][un] = 6
        if f[3] == ia or f[3] == 'N':
            d[na][un] = 6
        if f[3] == ic or f[3] == 'N':
            d[nc][un] = 6

    def orthorhombic_origin():
        """Fortran label 1399: orthorhombic groups.  Sets the
        origin-definition flags and returns the d-glide count."""
        nonlocal laueno
        laueno = 3
        # set up counts of the various types of mirrors
        im = ir = ia = ib = ic = id_ = i21 = 0
        if f2[0] == '2':
            # 2-axis along A; is it a 21?
            if f2[1] == '1':
                d[0][1] = 6
                d[0][2] = 6
                i21 = 4
        else:
            ir = 1
            if f2[0] == 'M':
                im = 4
            if f2[0] == 'B':
                ib = 1
            if f2[0] == 'C':
                ic = 1
            if f2[0] == 'D':
                id_ = 1
            if f3[0] == 'A' or f3[0] == 'N':
                d[0][0] = 6
            if f4[0] == 'A' or f4[0] == 'N':
                d[0][0] += 6
        if f3[0] == '2':
            # 2-axis along B; is it a 21?
            if f3[1] == '1':
                d[1][0] = 6
                d[1][2] = 6
                i21 += 2
        else:
            ir += 1
            if f3[0] == 'M':
                im += 2
            if f3[0] == 'A':
                ia = 1
            if f3[0] == 'C':
                ic += 1
            if f3[0] == 'D':
                id_ += 1
            if f2[0] == 'B' or f2[0] == 'N':
                d[1][1] = 6
            if f4[0] == 'B' or f4[0] == 'N':
                d[1][1] += 6
        if f4[0] == '2':
            # 2-axis along C; is it a 21?
            if f4[1] == '1':
                d[2][0] = 6
                d[2][1] = 6
                i21 += 1
        else:
            ir += 1
            if f4[0] == 'M':
                im += 1
            if f4[0] == 'A':
                ia += 1
            if f4[0] == 'B':
                ib += 1
            if f4[0] == 'D':
                id_ += 1
            if f2[0] == 'C' or f2[0] == 'N':
                d[2][2] = 6
            if f3[0] == 'C' or f3[0] == 'N':
                d[2][2] += 6
        if ir == 3:
            # 3 mirrors present; check for centering which may alter the
            # origin location
            if lcent == 2:
                # A-centered lattice
                if ib + ic == 1 and ia != 2:
                    d[1][1] += 6
                    d[2][2] += 6
            elif lcent == 3:
                # B-centered lattice
                if ia + ic == 1 and ib != 2:
                    d[0][0] += 6
                    d[2][2] += 6
            elif lcent == 4:
                # C-centered lattice
                if ia + ib == 1 and ic != 2:
                    d[0][0] += 6
                    d[1][1] += 6
            elif lcent == 5:
                # I-centered lattice: if only 1 glide plane shift mirrors by I
                if ia + ib + ic == 1:
                    d[0][0] += 6
                    d[1][1] += 6
                    d[2][2] += 6
        else:
            # less than 3 mirrors; set up the 2-axes locations
            if i21 in (4, 5, 7):
                d[0][1] = 0
            if i21 in (6, 7):
                d[0][2] = 0
            if i21 == 3:
                d[1][0] = 0
            if i21 in (2, 6, 7):
                d[1][2] = 0
            if i21 in (1, 3, 7):
                d[2][0] = 0
            if i21 == 5:
                d[2][1] = 0
            if im == 1 and i21 in (4, 2) and d[2][2] != 0:
                d[2][2] = 0
                d[2][1] += 6
            elif im == 2 and i21 in (4, 1) and d[1][1] != 0:
                d[1][1] = 0
                d[1][0] += 6
            elif im == 4 and i21 in (2, 1) and d[0][0] != 0:
                d[0][0] = 0
                d[0][2] += 6
        return id_

    # Now determine the Laue group and the unique axis if monoclinic
    if k == 2:
        # only 2 fields were read
        if f2[0] == '6':
            laueno = 11                     # 6/m
        elif f2[0] == '3':
            laueno = 8                      # 3bar
        elif f2[0] == '4':
            laueno = 4                      # 4/m
            if lcent >= 5:
                # I-centered or F-centered?
                if f2[3] != 'A' and f2[3] != 'D':
                    ier = tetragonal_tail(1553)
                else:
                    d[0][2] = 9
                    if lcent == 5:
                        d[1][2] = 3
            elif lcent == 4:
                # C-centered 4/m tetragonal
                if f2[2] == 'A' or f2[3] == 'A':
                    d[0][2] = 3
                    d[1][2] = 3
                    if f2[3] == 'A':
                        d[1][2] = 9
            elif f2[2] == 'N':
                d[0][2] = 6                 # P 4n/n * *
            elif f2[3] == 'N':
                d[1][2] = 6
        elif f2[0] == '1':
            laueno = 1                      # 1bar
        elif f2[0] == '5':
            ier = 5                         # bad 5-fold
        else:
            # 2/m, B-axis unique (short symbol)
            laueno = 2
            naxis = 2
            monoclinic_origin(1, 'A', 'C', 0, 1, 2, 1)
        return laueno, naxis, ier, 0, 0
    if k == 3:
        # only 3 fields were read; must be m3 cubic
        if f3[0] != '3':
            return 0, 0, 20, 0, 0
        laueno = 13
        if f2[1] == '1':
            d[1][0] = 6                     # B-axis flag if a 21 along A
        if f2[0] == 'B' or f2[0] == 'A':
            d[2][2] = 6                     # C-axis flag if an A-glide normal to C
        return laueno, 0, 0, 1, 0
    # four fields were read
    if f3[0] == '3':
        # m3m cubic
        laueno = 14
        if f2[0] == 'B' or f2[0] == 'A':
            d[2][2] = 6                     # C-axis translation if A or B normal to C
        if f2[0] == '4':
            # a 4n-axis specified
            if f2[1] == '-':
                # it is 4bar 3 *
                if f4[0] != 'M':
                    # it is not 4bar 3 m
                    if f4[0] == 'D':
                        # it is 4bar 3 d; I 4bar 3 d, we hope
                        if lcent != 5:
                            return laueno, 0, 21, 0, 0
                        d[0][2] = 9
                        d[1][2] = 3
                        d[2][2] = 9
                    else:
                        d[0][2] = d[1][2] = d[2][2] = 6
            elif f2[1] == '1':
                # 41-axis
                if lcent == 6:              # F 41 3 2
                    d[0][2] = 9
                    d[1][2] = 9
                    d[2][2] = 3
                else:                       # either P 41 3 2 or I 41 3 2
                    d[0][2] = 3
                    d[1][2] = 9
                    d[2][2] = 3
            elif f2[1] == '2':              # P 42 3 2
                d[0][2] = d[1][2] = d[2][2] = 6
            elif f2[1] == '3':              # 43 3 2
                d[0][2] = 9
                d[1][2] = 3
                d[2][2] = 9
        return laueno, 0, 0, 1, 0
    if f2[0] == '6':
        # hexagonal
        if f3[0] == '1' and f4[0] == '1':
            laueno = 11                     # 6/m
        else:
            laueno = 12                     # 6/mmm
    elif f2[0] == '3':
        # trigonal
        if f3[0] == '1':
            # P3**
            if f4[0] == '1':                # 31*
                laueno = 8                  # 3bar
            else:
                laueno = 10                 # 31m
        elif f4[0] != '1':
            laueno = 12                     # 6/mmm
        else:
            laueno = 9                      # 3m1
    elif f2[0] == '4':
        # tetragonal 4/mmm
        laueno = 5
        # if there is an N-glide normal to C place any mirror normal to A at 1/4
        if f2[2] == 'N' or f2[3] == 'N':
            d[0][0] = 6
        # if there is an A-glide normal to C place any mirror normal to (110) at 1/4
        if f2[2] == 'A' or f2[3] == 'A':
            d[1][1] = 3
        # if there is a 21 along B place it at x=1/4
        if f3[0] == '2' and f3[1] == '1':
            d[0][1] = 6
        # if there is a B- or N-glide normal to the A-axis shift the mirror by 1/4 along A
        if f3[0] == 'B' or f3[0] == 'N':
            d[0][0] += 6
        # if there is either a B- or N-glide normal to (110) shift the mirror by 1/4 along A
        if f4[0] == 'B' or f4[0] == 'N':
            d[1][1] += 3
        # if primitive and this is a 41, 42 or 43 and not 4n 21 2, set the
        # Z-location for 2-axes along (110)
        if lcent == 1 and f2[1] in ('1', '2', '3') and f3[1] != '1':
            d[2][0] = -3 * int(f2[1])
        # if fourth field is 21 and this is a 41, 42 or 43
        if f4[0] == '2' and f4[1] == '1' and f2[1] in ('1', '2', '3'):
            d[2][0] = 3 * int(f2[1])
        # set the Z-translation for 21-axes along B
        if f3[0] == '2' and f3[1] == '1' and f2[1] in ('1', '2', '3'):
            d[2][1] = 3 * int(f2[1])
        # place the D in F 4* D * at Y=7/8
        if _CODE[f3[0]] + _CODE[f2[2]] == 11 and lcent == 6:
            d[1][0] = 9
        # set M in F 4** * * at X=1/8 if a C along (110)
        if f4[0] == 'C' and lcent == 6:
            d[0][0] = 6
        if f2[1] == '-':
            lbl = 1556                      # this is a 4bar
        elif lcent > 1:
            lbl = 1553                      # the lattice is not primitive
        elif f2[2] == 'N' or f2[3] == 'N':
            lbl = 1552                      # N-glide normal to C
        elif f3[0] == '2' and f3[1] == '1':
            lbl = 1551                      # 21 along B
        elif f3[0] != 'N':
            lbl = 210
        elif f2[1] == '' or _CODE[f2[1]] > 15:
            lbl = 210
        else:
            lbl = 1551
        if lbl != 210:
            ier = tetragonal_tail(lbl)
    elif f2[0] == '1':
        if f3[0] == '1':
            if f4[0] == '1':
                laueno = 1                  # 1bar
            else:
                # C-axis unique monoclinic
                laueno = 2
                naxis = 3
                monoclinic_origin(3, 'A', 'B', 0, 2, 1, 2)
        elif f4[0] != '1':
            id_ = orthorhombic_origin()
        else:
            # B-axis unique monoclinic (full symbol)
            laueno = 2
            naxis = 2
            monoclinic_origin(2, 'A', 'C', 0, 1, 2, 1)
    elif f3[0] == '1':
        if f4[0] != '1':
            id_ = orthorhombic_origin()
        else:
            # A-axis unique monoclinic
            laueno = 2
            naxis = 1
            monoclinic_origin(1, 'B', 'C', 1, 0, 2, 0)
    else:
        # orthorhombic
        id_ = orthorhombic_origin()
    return laueno, naxis, ier, i209, id_


# --- the main interpreter (Fortran SGROUPNP) ----------------------------------

def _interpret(spg):
    """Generate a space group from its symbol.  Returns
    (laueno, naxis, ncent, lcent, npol, ops, ierr): the Laue group number,
    unique axis flag, centric flag, lattice centering flag, polar axis flag,
    the list of generated _Op symmetry operators and an error code."""
    laueno = naxis = ncent = lcent = npol = 0
    ops = []
    try:
        fields, k = _parse_symbol(spg)
        if k <= 1:
            raise _SgError(1)               # only 1 field was found
        if fields[0][0] not in _LATTICE:
            raise _SgError(2)               # first character was not P, A, B, C, F, I or R
        for fld in fields[1:]:
            if fld[0] == '-':
                # convert the -N notation to the Nbar notation (Fortran SGLPAK)
                if fld[1] not in ('1', '2', '3', '4', '5', '6'):
                    raise _SgError(4)
                fld[0], fld[1] = fld[1], '-'
        lcent = _LATTICE[fields[0][0]]
        d = [[0] * 3 for _ in range(3)]     # origin definition data
        i209 = 0                            # the body diagonal 3-axis flag
        id_ = 0                             # the number of d-glides
        if lcent != 7:
            laueno, naxis, ier, i209, id_ = _laue_group(fields, k, lcent, d)
            if ier:
                raise _SgError(ier)
        else:
            # rhombohedral lattice
            if fields[1][0] != '3':
                raise _SgError(3)           # make sure that there is a 3-axis
            if fields[k - 1][0] != 'R':
                # hexagonal axes, R centering; set laueno to 8 or 9
                if fields[k - 1][0] == 'H':
                    k -= 1
                laueno = k + 6
            else:
                # rhombohedral axes; delete R centering, set laueno to 6 or 7
                lcent = 1
                k -= 1
                laueno = k + 4
                i209 = 1
        cen = _CENTERING[lcent]             # list of lattice centering vectors
        ncv = len(cen)
        f2, f3, f4 = fields[1], fields[2], fields[3]

        def seed_mirror(i, mf, f):
            """Character was A, B, C, M or N: a mirror or axial glide plane."""
            t = [0, 0, 0]
            if mf == 1 and laueno <= 3:
                if k == 2:
                    # monoclinic B-axis unique: a B-axis mirror
                    if i == 2:
                        raise _SgError(9)
                    m = _M_B
                    t[1] = d[1][1]
                    if i in (1, 5):
                        t[0] = 6
                    if i in (3, 5):
                        t[2] = 6
                else:
                    # an A-axis mirror
                    if i == 1:
                        raise _SgError(8)
                    m = _M_A
                    t[0] = d[0][0]
                    if i in (2, 5):
                        t[1] = 6
                    if i in (3, 5):
                        t[2] = 6
            elif mf == 2 and laueno != 7:
                # third field and not a rhombohedral lattice
                if f2[0] in ('3', '6'):
                    # mirror normal to [100] in hex cell
                    m = _M_HEX
                    if i == 3:
                        t[2] = 6
                elif f2[0] == '4':
                    # it is not trigonal or hexagonal: an A-axis mirror
                    if i == 1:
                        raise _SgError(8)
                    m = _M_A
                    t[0] = d[0][0]
                    if i in (2, 5):
                        t[1] = 6
                    if i in (3, 5):
                        t[2] = 6
                else:
                    # a B-axis mirror
                    if i == 2:
                        raise _SgError(9)
                    m = _M_B
                    t[1] = d[1][1]
                    if i in (1, 5):
                        t[0] = 6
                    if i in (3, 5):
                        t[2] = 6
            elif mf == 3 or laueno > 3:
                if (mf == 3 or laueno == 7) and \
                        (f3[0] == '3' or f2[0] == '4' or f2[0] in ('3', '6')):
                    # it is not cubic or tetragonal:
                    # a diagonal mirror normal to [-110]
                    m = _M_DIAG
                    t[0] = d[1][1]
                    t[1] = -d[1][1]
                    if i in (3, 5):
                        t[2] = 6
                    if (laueno == 7 and i == 3) or i < 3 or i > 4:
                        if lcent in (6, 4):
                            # F- or C-centered tetragonal: glides are 1/4,1/4
                            t[0] += 3
                            t[1] += 3
                        else:
                            t[0] += 6
                            t[1] += 6
                else:
                    # a C-axis mirror
                    if i == 3:
                        raise _SgError(10)
                    m = _M_C
                    t[2] = d[2][2]
                    if i in (1, 5):
                        t[0] = 6
                    if i in (2, 5):
                        t[1] = 6
                    if mf == 1 and f2[0] == '6' and f2[1] == '3':
                        # if this is a 63-axis the mirror is at 1/4
                        t[2] = 6
            return _Op(m, t)

        def seed_d_glide(mf):
            """A d-glide type mirror."""
            if lcent <= 1:
                raise _SgError(11)
            t = [0, 0, 0]
            icv = 1
            if mf == 1 and laueno <= 3:
                if k == 2:
                    if ncv == 4:
                        icv = 2
                    m = _M_B
                    t[0] = cen[icv][0] // 2
                    t[2] = cen[icv][2] // 2
                else:
                    m = _M_A
                    if id_ == 2:
                        t[0] = 3
                    t[1] = cen[icv][1] // 2
                    t[2] = cen[icv][2] // 2
            elif mf == 2:
                if ncv == 4:
                    icv = 2
                m = _M_B
                t[0] = cen[icv][0] // 2
                if id_ == 2:
                    t[1] = 3
                if laueno == 5:
                    t[1] = d[1][0]
                t[2] = cen[icv][2] // 2
            elif mf == 3 or laueno > 3:
                if mf == 3 and (f2[0] == '4' or f3[0] == '3'):
                    # cubic or tetragonal: d-glide along diagonal
                    m = _M_DIAG
                    if f3[0] == '2':
                        t = [0, 6, 3]
                    else:
                        t = [3, 3, 3]
                else:
                    if ncv == 4:
                        icv = 3
                    m = _M_C
                    t[0] = cen[icv][0] // 2
                    t[1] = cen[icv][1] // 2
                    if id_ == 2:
                        t[2] = 3
            return _Op(m, t)

        def seed_two_fold(mf, f):
            """A 2-fold rotation axis."""
            if f[1] == '-':
                raise _SgError(19)          # we will not allow a -2 axis
            t = [0, 0, 0]
            if mf == 1:
                # first rotation operator
                if k == 2:
                    # rotation about the B-axis
                    m = _2_B
                    t[0] = d[0][1]
                    t[2] = d[2][1]
                    if f[1] == '1':
                        t[1] = 6
                else:
                    # rotation about the A-axis
                    m = _2_A
                    t[1] = d[1][0]
                    t[2] = d[2][0]
                    if f[1] in ('1', '3'):
                        t[0] = 6
            elif mf == 2:
                # second rotation operator
                if laueno == 7:
                    m = _2_1M10             # 2-axis along [1-10]
                elif f2[0] == '6' and f4[0] != '1':
                    m = _2_100H             # 2-axis along [100] (P 6n22 groups)
                elif f2[0] == '3':
                    # 2-axis along [110] trig; also used for the P 3n21 groups
                    m = _2_110
                    t[0] = d[1][0]
                    if f[1] == '1':
                        t[0] += 6
                    t[1] = -d[1][0]
                    t[2] = d[2][0]
                else:
                    # not a hexagonal or trigonal space group:
                    # rotation about the B-axis
                    m = _2_B
                    if f2[0] == 'M' and f4[0] == 'N':
                        t[0] = 6
                    else:
                        t[0] = d[0][1]
                    t[2] = d[2][1]
                    if f[1] == '1':
                        t[1] = 6
            elif mf == 3:
                if _CODE[f2[0]] >= 14 or f3[0] == '3':
                    if f2[0] == '4':
                        # 2-axis along [110] tetrag
                        m = _2_110
                        t[0] = d[1][0]
                        if f[1] == '1':
                            t[0] += 6
                        t[1] = -d[1][0]
                        t[2] = d[2][0]
                    else:
                        m = _2_210          # 2-axis along [210]
                else:
                    # 2-fold rotation about the C-axis
                    m = _2_C
                    t[0] = d[0][2]
                    t[1] = d[1][2]
                    if f[1] in ('1', '3'):
                        t[2] = 6
                    if f[1] == '5':
                        t[2] = 6
            return _Op(m, t)

        ops.append(_Op(_IDENTITY))
        slot = _Op(_ZERO)                   # stale content of the next slot
        if i209:
            # cubic or rhombohedral cell: generate z,x,y and y,z,x
            ops.append(_Op(_ZXY))
            ops.append(_Op(_YZX))
        ierx = 0
        for mf in range(1, k):
            # decode the last 3 fields of the symbol
            f = fields[mf]
            ifld = 0
            while ifld < 4 and f[ifld]:
                if ifld > 0:
                    # skip to the token after a '/' (for 2/m style fields)
                    while ifld <= 2 and f[ifld] != '/':
                        if not f[ifld]:
                            ifld = 3
                        else:
                            if f[ifld] in _LETTERS:
                                raise _SgError(16)
                            ifld += 1
                    ifld += 1
                    if ifld < 4 and not f[ifld]:
                        raise _SgError(17)
                if ifld < 4:
                    i = _OP.get(f[ifld])
                    if i is None:
                        raise _SgError(7)
                    ndelt = True
                    nxi = len(ops)          # set first matrix pointer
                    op = None
                    if i <= 5:
                        op = seed_mirror(i, mf, f)
                    elif i == 6:
                        op = seed_d_glide(mf)
                    elif i == 7:
                        # 1-fold axis
                        ndelt = False
                        if f[1] == '-':
                            ncent = 1       # we have a center of symmetry
                            ifld += 1
                    elif i == 8:
                        op = seed_two_fold(mf, f)
                    elif i == 9:
                        # 3-fold axis
                        if mf == 1 and laueno > 7:
                            t = [0, 0, 0]
                            if f[1] == '1':
                                t[2] = 4
                            if f[1] == '2':
                                t[2] = 8
                            if f[1] == '-':
                                ncent = 1
                                ifld += 1
                            op = _Op(_3FOLD, t)
                        elif mf == 2 or laueno <= 7:
                            ndelt = False
                            if f[1] == '-':
                                ncent = 1
                                ifld += 1
                        else:
                            raise _SgError(25)
                    elif i == 10:
                        # four fold axis
                        if mf != 1:
                            raise _SgError(12)
                        if f[1] == '-':
                            # 4-bar axis
                            op = _Op(_4BAR, [d[0][2], d[1][2], d[2][2]])
                            ifld += 1
                        else:
                            # 4-axis
                            t = [d[0][2], d[1][2], 0]
                            if f[1] == '1':
                                t[2] = 3    # 41 axis
                            if f[1] == '2':
                                t[2] = 6    # 42 axis
                            if f[1] == '3':
                                t[2] = 9    # 43 axis
                            op = _Op(_4FOLD, t)
                    elif i == 12:
                        # 6-axis
                        if mf != 1:
                            raise _SgError(13)
                        if f[1] == '-':
                            # 6-bar operation
                            t = [0, 0, 0]
                            if f3[0] == 'C' or f4[0] == 'C':
                                t[2] = 6
                            op = _Op(_6BAR, t)
                            ifld += 1
                        else:
                            # 6 operation
                            t = [0, 0, 0]
                            if f[1] in ('1', '2', '3', '4', '5'):
                                t[2] = 2 * int(f[1])
                            op = _Op(_6FOLD, t)
                    # i in (11, 13, 14, 15): no operator is generated
                    if ndelt:
                        if op is not None:
                            slot = op
                        added, cent, ier = _try_add(ops, slot, mf, lcent, laueno)
                        if cent:
                            ncent = 1
                        if ier:
                            ierx = ier
                        if added:
                            if len(ops) >= 25:
                                raise _SgError(14)  # should never be more than 24
                            nxl = len(ops) - 1
                            # repeat this loop until no new matrices
                            while nxi <= nxl:
                                for nx in range(nxi, nxl + 1):
                                    for m1 in range(1, nx + 1):
                                        # apply NX to M1 to generate a matrix
                                        slot = _compose(ops[m1], ops[nx])
                                        added, cent, ier = _try_add(
                                            ops, slot, mf, lcent, laueno)
                                        if cent:
                                            ncent = 1
                                        if ier:
                                            ierx = ier
                                        if added and len(ops) >= 25:
                                            raise _SgError(14)
                                nxi = nxl + 1       # first matrix to first new matrix
                                nxl = len(ops) - 1  # last matrix
                ifld += 1
        # determine presence of indeterminate origin
        npx, npy, npz, npxyz, npyxz = 1, 2, 4, 0, 1
        for op in ops:
            if op.m[0, 0] <= 0:
                npx = 0                     # origin is defined along X
            if op.m[1, 1] <= 0:
                npy = 0                     # origin is defined along Y
            if op.m[2, 2] <= 0:
                npz = 0                     # origin is defined along Z
            if op.m[0, 2] > 0:
                npxyz = 8                   # there is a 3-axis along [1,1,1]
            if op.m[0, 2] < 0:
                npyxz = 0                   # origin is defined along [1,1,1]
        # set the indeterminate origin flag
        npol = (npx + npy + npz + npxyz * npyxz) * (1 - ncent)
        if lcent == 7:
            naxis = -1
        if ierx:
            return laueno, 4, ncent, lcent, npol, ops, ierx
        return laueno, naxis, ncent, lcent, npol, ops, 0
    except _SgError as err:
        return laueno, 4, ncent, lcent, 0, [], err.args[0]


# --- the public entry point (Fortran SGFORPY in pyspg.for) ---------------------

def sgforpy(spg):
    """Space group access routine: generate the symmetry operators for a
    space-group symbol.

    :param str spg: the space-group symbol (e.g. ``'P 21 21 21'``); at most
        20 characters are used, lower case is accepted and a tab is treated
        as a space.
    :returns: the list ``[LAUE, SGINV, SGLATT, SGUNIQ, SGPOL, SGNOPS,
        SGMTRX, SGTRNS, SGGEN, IERR]`` with entries:

        * **LAUE** (*int*) -- Laue group number: 1=-1, 2=2/m, 3=mmm, 4=4/m,
          5=4/mmm, 6=R3R, 7=R3MR, 8=3, 9=3m1, 10=31m, 11=6/m, 12=6/mmm,
          13=m3, 14=m3m
        * **SGINV** (*int*) -- centrosymmetry flag: 0 = acentric,
          1 = centric
        * **SGLATT** (*int*) -- lattice centering: 1=P, 2=A, 3=B, 4=C, 5=I,
          6=F, 7=R
        * **SGUNIQ** (*int*) -- unique axis for monoclinic groups: 1=a, 2=b,
          3=c; -1 for rhombohedral in the hexagonal setting; 4 on error exit
        * **SGPOL** (*int*) -- polar axis flag: 1=x, 2=y, 3=xy, 4=z, 5=xz,
          6=yz, 7=xyz, 8=111
        * **SGNOPS** (*int*) -- number of symmetry operators generated
        * **SGMTRX** (*numpy.ndarray*) -- (SGNOPS,3,3) rotation matrices
        * **SGTRNS** (*numpy.ndarray*) -- (SGNOPS,3) translation vectors
          (fractional)
        * **SGGEN** (*numpy.ndarray*) -- (SGNOPS,) symmetry operation flags
        * **IERR** (*int*) -- error flag: 0 = no error
    """
    laueno, naxis, ncent, lcent, npol, ops, ierr = _interpret(spg)
    nops = len(ops)
    sgmtrx = np.zeros((nops, 3, 3))
    sgtrns = np.zeros((nops, 3))
    sggen = np.zeros(nops, dtype=int)
    for k, op in enumerate(ops):
        sgmtrx[k] = op.m
        sgtrns[k] = op.t / 12.0
        sggen[k] = _OPFLAG.get(op.key, -1)
    return [laueno, ncent, lcent, naxis, npol, nops, sgmtrx, sgtrns, sggen,
            ierr]


if __name__ == '__main__':
    import sys
    for line in sys.stdin:
        spg = line.strip()
        if not spg:
            continue
        (laue, sginv, sglatt, sguniq, sgpol, sgnops, sgmtrx, sgtrns, sggen,
         ierr) = sgforpy(spg)
        print(spg)
        if ierr:
            print('ERROR with', spg, 'IERR =', ierr)
            continue
        print(laue, sginv, sglatt, sguniq, 0)
        print(sgnops)
        for k in range(sgnops):
            for i in range(3):
                print(list(sgmtrx[k][i]), '+', sgtrns[k][i])
            print(sggen[k])
