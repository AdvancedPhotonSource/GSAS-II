#!/usr/bin/env python
##############################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Chris Farros, Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""PdfFit class for fitting pdf data to a model."""

from __future__ import print_function

import itertools
import numbers

import six

# Load pdffit2 and output modules to the current namespace.
# Note that "import diffpy.pdffit2.output as output" would
# crash with AttributeError when executed during imports of
# parent packages.
try:
    from diffpy.pdffit2 import pdffit2
except:
    import pdffit2     # added for GSAS-II to relocate binary file
from diffpy.pdffit2 import output

# helper routines

def _format_value_std(value, stdev):
    """Convert value to a string with standard deviation in brackets.

    value -- a number
    stdev -- standard deviation.  Ignored when small compared to value.

    Return string.
    """
    if stdev > abs(value)*1e-8:
        s = "%g (%g)" % (value, stdev)
    elif str(stdev) == 'nan':
        s = "%g (NaN)" % value
    else:
        s = "%g" % value
    return s


def _format_bond_length(dij, ddij, ij1, symij):
    """Return string with formatted bond length info for a pair of atoms.

    dij     -- distance between atoms i and j
    ddij    -- standard deviation of dij.  Ignored when small relative to dij.
    ij1     -- tuple of atom indices starting at 1
    symij   -- tuple of atom symbols

    Return formatted string.
    """
    w_smbidx = 10
    w_equals = 30
    s0 = "%s (#%i)" % (symij[0], ij1[0])
    s1 = "%s (#%i)" % (symij[1], ij1[1])
    leader0 = "   " + s0.ljust(w_smbidx) + " -   " + s1 + " "
    leader1 = leader0.ljust(w_equals) + '=   '
    s = leader1 + _format_value_std(dij, ddij) + " A"
    return s


def _convertCallable(var):
    """Convert an object to the result of its call when callable.

    var -- string or callable object that returns string

    Return var or var().
    """
    if callable(var):
        rv = var()
    else:
        rv = var
    return rv


# constants

__intro_message__ = """
******************************************************************************
*                      P D F F I T   Version   %(version)s                   *
*                                              %(date)s                      *
* -------------------------------------------------------------------------- *
* (c) 1998-2007 Trustees of the Michigan State University.                   *
* (c) 2008-%(year)s Trustees of the Columbia University                      *
*               in the city of New York.                                     *
*                                                                            *
* Authors:                                                                   *
*     Thomas Proffen        -   Email: tproffen@lanl.gov                     *
*     Jacques Bloch         -   Email: bloch@pa.msu.edu                      *
*     Christopher Farrow    -   Email: clf2121@columbia.edu                  *
*     Pavol Juhas           -   Email: pjuhas@bnl.gov                        *
*     Simon Billinge        -   Email: sb2896@columbia.edu                   *
******************************************************************************
"""


##############################################################################


class PdfFit(object):
    """Create PdfFit object."""

    # constants and enumerators from pdffit.h:
    # selection of all atoms
    selalias = { 'ALL' : -1 }
    # constraint type identifiers
    FCON = { 'USER' : 0, 'IDENT' : 1, 'FCOMP' : 2, 'FSQR' : 3 }
    # scattering type identifiers
    Sctp = { 'X' : 0, 'N' : 1 }

    def _exportAll(self, namespace):
        """ _exportAll(self, namespace) --> Export all 'public' class methods
            into namespace.

        This function allows for a module-level PdfFit object which doesn't have
        to be referenced when calling a method. This function makes old (python)
        scripts compatible with this class. At the top of the script, create a
        pdffit object, and then call this method. Usually, namespace = locals().
        """
        # string aliases (var = "var")
        for a in itertools.chain(self.selalias, self.FCON, self.Sctp):
            exec("%s = %r" % (a, a), namespace)
        public = [ a for a in dir(self) if "__" not in a and a not in
                ["_handle", "_exportAll", "selalias", "FCON", "Sctp" ] ]
        for funcname in public:
            namespace[funcname] = getattr(self, funcname)
        return

    def intro():
        """Show introductory message.
        """
        import re
        from diffpy.pdffit2 import __version__, __date__
        date = __date__[:10]
        d = {'version' : __version__,  'date' : date,
             'year' : date[:4] or '2019'}
        msg = __intro_message__ % d
        filler = lambda mx : (mx.group(0).rstrip(' *').ljust(77) + '*')
        msg_ljust = re.sub('(?m)^(.{1,77}|.{79}.*)$', filler, msg)
        print(msg_ljust, file=output.stdout)
        return
    intro = staticmethod(intro)


    def add_structure(self, stru):
        """add_structure(stru) --> Add new structure to PdfFit instance.

        stru -- instance of Structure class from diffpy.structure.

        No return value.
        Raises pdffit2.structureError when stru contains unknown
        atom species.
        """
        s = stru.writeStr('pdffit')
        self.read_struct_string(s)
        return


    def read_struct(self, struct):
        """read_struct(struct) --> Read structure from file into memory.

        struct  -- name of file from which to read structure

        Raises:
            pdffit2.calculationError when a lattice cannot be created from the
            given structure
            pdffit2.structureError when a structure file is malformed
            IOError when the file cannot be read from the disk
        """
        pdffit2.read_struct(self._handle, struct)
        self.stru_files.append(struct)
        return


    def read_struct_string(self, struct, name=""):
        """read_struct_string(struct, name = "") --> Read structure from
        a string into memory.

        struct  -- string containing the contents of the structure file
        name    -- tag with which to label structure

        Raises:
            pdffit2.calculationError when a lattice cannot be created from the
            given structure
            pdffit2.structureError when a structure file is malformed
        """
        pdffit2.read_struct_string(self._handle, struct)
        self.stru_files.append(name)
        return


    def read_data(self, data, stype, qmax, qdamp):
        """read_data(data, stype, qmax, qdamp) --> Read pdf data from file into
        memory.

        data    -- name of file from which to read data
        stype   -- 'X' (xray) or 'N' (neutron)
        qmax    -- Q-value cutoff used in PDF calculation.
                   Use qmax=0 to neglect termination ripples.
        qdamp   -- instrumental Q-resolution factor

        Raises: IOError when the file cannot be read from disk
        """
        pdffit2.read_data(self._handle, data, six.b(stype), qmax, qdamp)
        self.data_files.append(data)
        return


    def read_data_string(self, data, stype, qmax, qdamp, name=""):
        """read_data_string(data, stype, qmax, qdamp, name = "") --> Read
        pdf data from a string into memory.

        data    -- string containing the contents of the data file
        stype   -- 'X' (xray) or 'N' (neutron)
        qmax    -- Q-value cutoff used in PDF calculation.
                   Use qmax=0 to neglect termination ripples.
        qdamp   -- instrumental Q-resolution factor
        name    -- tag with which to label data
        """
        pdffit2.read_data_string(self._handle, data, six.b(stype), qmax,
                qdamp, name)
        name = data
        self.data_files.append(name)
        return


    def read_data_lists(self, stype, qmax, qdamp, r_data, Gr_data,
            dGr_data = None, name = "list"):
        """read_data_lists(stype, qmax, qdamp, r_data, Gr_data, dGr_data =
        None, name = "list") --> Read pdf data into memory from lists.

        All lists must be of the same length.
        stype       -- 'X' (xray) or 'N' (neutron)
        qmax        -- Q-value cutoff used in PDF calculation.
                       Use qmax=0 to neglect termination ripples.
        qdamp       -- instrumental Q-resolution factor
        r_data      -- list of r-values
        Gr_data     -- list of G(r) values
        dGr_data    -- list of G(r) uncertainty values
        name        -- tag with which to label data

        Raises: ValueError when the data lists are of different length
        """
        pdffit2.read_data_arrays(self._handle, six.b(stype), qmax, qdamp,
                r_data, Gr_data, dGr_data, name)
        self.data_files.append(name)
        return


    def pdfrange(self, iset, rmin, rmax):
        """pdfrange(iset, rmin, rmax) --> Set the range of the fit.

        iset    -- data set to consider
        rmin    -- minimum r-value of fit
        rmax    -- maximum r-value of fit

        Raises: ValueError for bad input values
        """
        pdffit2.pdfrange(self._handle, iset, rmin, rmax)
        return


    def reset(self):
        """reset() --> Clear all stored fit, structure, and parameter data."""
        self.stru_files = []
        self.data_files = []
        pdffit2.reset(self._handle);
        return


    def alloc(self, stype, qmax, qdamp, rmin, rmax, bin):
        """alloc(stype, qmax, qdamp, rmin, rmax, bin) --> Allocate space
        for a PDF calculation.

        The structure from which to calculate the PDF must first be imported
        with the read_struct() or read_struct_string() method.
        stype   -- 'X' (xray) or 'N' (neutron)
        qmax    -- Q-value cutoff used in PDF calculation.
                   Use qmax=0 to neglect termination ripples.
        qdamp   -- instrumental Q-resolution factor
        rmin    -- minimum r-value of calculation
        rmax    -- maximum r-value of calculation
        bin     -- number of data points in calculation

        Raises:
            ValueError for bad input values
            pdffit.unassignedError when no structure has been loaded
        """
        pdffit2.alloc(self._handle, six.b(stype), qmax, qdamp, rmin,
                rmax, bin)
        return


    def calc(self):
        """calc() --> Calculate the PDF of the imported structure.

        Space for the calculation must first be allocated with the alloc()
        method.

        Raises:
            pdffit2.calculationError when allocated space cannot
            accomodate calculation
            pdffit.unassignedError when space for calculation has not been
            allocated
        """
        pdffit2.calc(self._handle)
        return


    def refine(self, toler=0.00000001):
        """refine(toler = 0.00000001) --> Fit the theory to the imported data.

        toler   --  tolerance of the fit

        Raises:
            pdffit2.calculationError when the model pdf cannot be calculated
            pdffit2.constraintError when refinement fails due to bad
            constraint
            pdffit2.unassigedError when a constraint used but never initialized
            using setpar()
        """
        step = 0
        finished = 0
        while not finished:
            finished = pdffit2.refine_step(self._handle, toler)
            step += 1
        return


    def refine_step(self, toler=0.00000001):
        """refine_step(toler = 0.00000001) --> Run a single step of the fit.

        toler   --  tolerance of the fit

        Raises:
            pdffit2.calculationError when the model pdf cannot be calculated
            pdffit2.constraintError when refinement fails due to bad
            constraint
            pdffit2.unassigedError when a constraint used but never initialized
            using setpar()

        Returns: 1 (0) if refinement is (is not) finished
        """
        self.finished = pdffit2.refine_step(self._handle, toler)
        return self.finished


    def save_pdf(self, iset, fname):
        """save_pdf(iset, fname) --> Save calculated or fitted PDF to file.

        iset    -- data set to save

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if the data set is undefined
        """
        pdffit2.save_pdf(self._handle, iset, fname)
        return


    def save_pdf_string(self, iset):
        """save_pdf_string(iset) --> Save calculated or fitted PDF to string.

        iset    -- data set to save

        Raises:
            pdffit2.unassignedError if the data set is undefined

        Returns: string containing contents of save file
        """
        pdffilestring = pdffit2.save_pdf(self._handle, iset, "")
        return pdffilestring


    def save_dif(self, iset, fname):
        """save_dif(iset, fname) --> Save data and fitted PDF difference to
        file.

        iset    -- data set to save

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if the data set is undefined
        """
        pdffit2.save_dif(self._handle, iset, fname)
        return


    def save_dif_string(self, iset):
        """save_dif_string(iset) --> Save data and fitted PDF difference to
        string.

        iset    -- data set to save

        Raises:
            pdffit2.unassignedError if the data set is undefined

        Returns: string containing contents of save file
        """
        diffilestring = pdffit2.save_dif(self._handle, iset, "")
        return diffilestring


    def save_res(self, fname):
        """save_res(fname) --> Save fit-specific data to file.

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if there is no refinement data to save
        """
        pdffit2.save_res(self._handle, fname)
        return


    def save_res_string(self):
        """save_res_string() --> Save fit-specific data to a string.

        Raises:
            pdffit2.unassignedError if there is no refinement data to save

        Returns: string containing contents of save file
        """
        resfilestring = pdffit2.save_res(self._handle, "")
        return resfilestring


    def get_structure(self, ip):
        """get_structure(ip) --> Get a copy of specified phase data.

        ip -- index of existing PdfFit phase starting from 1

        Return Structure object from diffpy.structure.
        Raise pdffit2.unassignedError if phase ip is undefined.
        """
        from diffpy.structure import PDFFitStructure
        s = self.save_struct_string(ip)
        stru = PDFFitStructure()
        stru.readStr(s, 'pdffit')
        return stru


    def save_struct(self, ip, fname):
        """save_struct(ip, fname) --> Save structure resulting from fit
        to file.

        ip    -- phase to save

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if the data set is undefined
        """
        pdffit2.save_struct(self._handle, ip, fname)
        return


    def save_struct_string(self, ip):
        """save_struct(ip) --> Save structure resulting from fit to string.

        ip    -- phase to save

        Raises:
            pdffit2.unassignedError if phase ip is undefined.

        Returns: string containing contents of save file
        """
        structfilestring = pdffit2.save_struct(self._handle, ip, "")
        return structfilestring


    def show_struct(self, ip):
        """show_struct(ip) --> Print structure resulting from fit.

        ip    -- phase to display

        Raises: pdffit2.unassignedError if the phase is undefined
        """
        pdffit2.show_struct(self._handle, ip)
        return


    def constrain(self, var, par, fcon=None):
        """constrain(var, par[, fcon]) --> Constrain a variable to a parameter.

        A variable can be constrained to a number or equation string.
        var     -- variable to constrain, such as x(1)
        par     -- parameter which to constrain the variable. This can be
                   an integer or an equation string containing a reference
                   to another parameter. Equation strings use standard c++
                   syntax. The value of a constrained parameter is accessed
                   as @p in an equation string, where p is the parameter.
                   e.g.
                   >>>  constrain(x(1), 1)
                   >>>  constrain(x(2), "0.5+@1")
        fcon    -- 'USER', 'IDENT', 'FCOMP', or 'FSQR'
                   this is an optional parameter, and I don't know how it is
                   used!

        Raises:
            pdffit2.constraintError if a constraint is bad
            pdffit2.unassignedError if variable does not yet exist
            ValueError if variable index does not exist (e.g. lat(7))
        """
        var_ref = self.__getRef(var)
        varnc = _convertCallable(var)
        if fcon:
            fc = self.FCON[fcon]
            pdffit2.constrain_int(self._handle, var_ref, varnc, par, fc)
        elif isinstance(par, six.string_types):
            pdffit2.constrain_str(self._handle, var_ref, varnc, par)
        else:
            pdffit2.constrain_int(self._handle, var_ref, varnc, par)
        return


    def setpar(self, par, val):
        """setpar(par, val) --> Set value of constrained parameter.

        val     --  Either a numerical value or a reference to a variable

        Raises:
            pdffit2.unassignedError when variable is yet to be assigned
        """
        # people do not use parenthesis, e.g., "setpar(3, qdamp)"
        # in such case val is a reference to PdfFit method
        val = _convertCallable(val)
        try:
            val = float(val)
            pdffit2.setpar_dbl(self._handle, par, val)
        except ValueError:
            var_ref = self.__getRef(val)
            pdffit2.setpar_RV(self._handle, par, var_ref)
        return


    def setvar(self, var, val):
        """setvar(var, val) --> Set the value of a variable.

        Raises:
            pdffit2.unassignedError if variable does not yet exist
            ValueError if variable index does not exist (e.g. lat(7))
        """
        var_ref = self.__getRef(var)
        pdffit2.setvar(self._handle, var_ref, val)
        return


    def getvar(self, var):
        """getvar(var) --> Get stored value of a variable.

        Raises:
            pdffit2.unassignedError if variable does not yet exist
            ValueError if variable index does not exist (e.g. lat(7))
        """
        var_ref = self.__getRef(var)
        retval = pdffit2.getvar(self._handle, var_ref)
        return retval


    def getrw(self):
        """getrw() --> Get normalized total error of the fit rw.

        getrw calculates total fit error summed for all datasets in the fit.

        Return float.
        """
        rw = pdffit2.getrw(self._handle)
        return rw


    def getcrw(self):
        """getcrw() --> Get cumulative Rw for the current dataset.

        Cumulative Rw is a list of Rw partial sums cost values evaluated against
        observed PDF data in the error sums evaluated against
        the r-points in the fit.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of crw points, equidistant in r or empty list
        if the refine function has not been called yet.
        """
        crw = pdffit2.getcrw(self._handle)
        return crw


    def getR(self):
        """getR() --> Get r-points used in the fit.

        This function should only be called after data has been loaded or
        calculated. Before a refinement, the list of r-points will reflect the
        data. Afterwords, they will reflect the fit range.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of equidistance r-points used in fit.
        """
        R = pdffit2.getR(self._handle)
        return R


    def getpdf_fit(self):
        """getpdf_fit() --> Get fitted PDF.

        This function should only be called after a refinement or refinement
        step has been done.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of fitted points, equidistant in r.
        """
        pdfdata = pdffit2.getpdf_fit(self._handle)
        return pdfdata


    def getpdf_obs(self):
        """getpdf_obs() --> Get observed PDF.

        This function should only be called after data has been loaded or
        calculated. Before a refinement, the list of r-points will reflect the
        data. Afterwords, they will reflect the fit range.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of data points, equidistant in r.
        """
        pdfdata = pdffit2.getpdf_obs(self._handle)
        return pdfdata


    def getpdf_diff(self):
        """Obtain difference between observed and fitted PDF.

        This function should only be called after data has been loaded or
        calculated. Before a refinement, the list of r-points will reflect the
        data. Afterwords, they will reflect the fit range.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of data points, equidistant in r.
        """
        Gdiff = pdffit2.getpdf_diff(self._handle)
        return Gdiff


    def get_atoms(self, ip=None):
        """get_atoms() --> Get element symbols of all atoms in the structure.

        ip -- index of phase to get the elements from (starting from 1)
              when ip is not given, use current phase

        This function should only be called after a structure has been loaded.

        Raises: pdffit2.unassignedError if no structure exists

        Returns: List of atom names in structure.
        """
        if ip is None:  rv = pdffit2.get_atoms(self._handle)
        else:           rv = pdffit2.get_atoms(self._handle, ip)
        return rv


    def get_atom_types(self, ip=None):
        """get_atom_types() --> Ordered unique element symbols in the structure.

        ip -- index of phase to get the elements from (starting from 1)
              when ip is not given, use current phase

        This function should only be called after a structure has been loaded.

        Raises:
            pdffit2.unassignedError if no structure exists

        Returns: List of unique atom symbols as they occur in structure.
        """
        if ip is None:  rv = pdffit2.get_atom_types(self._handle)
        else:           rv = pdffit2.get_atom_types(self._handle, ip)
        return rv


    def getpar(self, par):
        """getpar(par) --> Get value of parameter.

        Raises: ValueError if parameter does not exists
        """
        return pdffit2.getpar(self._handle, par)


    def fixpar(self, par):
        """fixpar(par) --> Fix a parameter.

        Fixed parameters are not fitted in a refinement. Passed parameter
        can be 'ALL', in which case all parameters are fixed.

        Raises: pdffit.unassignedError when parameter has not been assigned
        """
        if isinstance(par, six.string_types) and par.upper() in self.selalias:
            par = self.selalias[par.upper()]
        pdffit2.fixpar(self._handle, par)
        return


    def freepar(self, par):
        """freepar(par) --> Free a parameter.

        Freed parameters are fitted in a refinement. Passed parameter
        can be 'ALL', in which case all parameters are freed.

        Raises: pdffit.unassignedError when parameter has not been assigned
        """
        if isinstance(par, six.string_types) and par.upper() in self.selalias:
            par = self.selalias[par.upper()]
        pdffit2.freepar(self._handle, par)
        return


    def setphase(self, ip):
        """setphase(ip) --> Switch to phase ip.

        ip  -- index of the phase starting at 1.

        All parameters assigned after this method is called refer only to the
        current phase.

        Raises: pdffit.unassignedError when phase does not exist
        """
        pdffit2.setphase(self._handle, ip)
        return


    def setdata(self, iset):
        """setdata(iset) --> Set the data set in focus.

        iset -- integer index of data set starting at 1.

        Raises: pdffit.unassignedError when data set does not exist
        """
        pdffit2.setdata(self._handle, iset)
        return


    def psel(self, ip):
        """psel(ip) --> Include phase ip in calculation of total PDF

        psel('ALL')     selects all phases for PDF calculation.

        Raises: pdffit2.unassignedError if selected phase does not exist
        """
        if isinstance(ip, six.string_types) and ip.upper() in self.selalias:
            ip = self.selalias[ip.upper()]
        pdffit2.psel(self._handle, ip)
        return


    def pdesel(self, ip):
        """pdesel(ip) --> Exclude phase ip from calculation of total PDF.

        pdesel('ALL')   excludes all phases from PDF calculation.

        Raises: pdffit2.unassignedError if selected phase does not exist
        """
        if isinstance(ip, six.string_types) and ip.upper() in self.selalias:
            ip = self.selalias[ip.upper()]
        pdffit2.pdesel(self._handle, ip)
        return


    def selectAtomType(self, ip, ijchar, symbol, flag):
        """Configure partial PDF - mark the specified atom type in phase ip
        as included or excluded as a first or second in pair for distance
        evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair
        symbol  -- element symbol
        flag    -- bool flag, True for selection, False for exclusion

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError for invalid value of ijchar
        """
        pdffit2.selectAtomType(self._handle, ip, six.b(ijchar), symbol, flag)
        return


    def selectAtomIndex(self, ip, ijchar, aidx, flag):
        """Configure partial PDF - mark the atom of given index in phase ip
        as included or excluded as a first or second in pair for distance
        evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair
        aidx    -- integer index of atom starting at 1
        flag    -- bool flag, True for selection, False for exclusion

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if atom index or ijchar are invalid
        """
        pdffit2.selectAtomIndex(self._handle, ip, six.b(ijchar), aidx, flag)
        return


    def selectAll(self, ip, ijchar):
        """Configure partial PDF - include all atoms of phase ip as first or
        second element in pair for distance evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if ijchar is invalid
        """
        pdffit2.selectAll(self._handle, ip, six.b(ijchar))
        return


    def selectNone(self, ip, ijchar):
        """Configure partial PDF - exclude all atoms of phase ip from first
        or second element of pair distance evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if ijchar is invalid
        """
        pdffit2.selectNone(self._handle, ip, six.b(ijchar))
        return


    def bang(self, i, j, k):
        """bang(i, j, k) --> Show bond angle defined by atoms i, j, k.

        No return value.  Use bond_angle() to get the result.

        Raises: ValueError if selected atom(s) does not exist
                pdffit.unassignedError when no structure has been loaded
        """
        angle, stdev = pdffit2.bond_angle(self._handle, i, j, k)
        # indices should be already checked here by bond_angle
        atom_symbols = self.get_atoms()
        leader = "   %s (#%i) - %s (#%i) - %s (#%i)   =   " % \
                (atom_symbols[i-1], i, atom_symbols[j-1], j,
                 atom_symbols[k-1], k)
        s = leader + _format_value_std(angle, stdev) + " degrees"
        print(s, file=output.stdout)
        return


    def bond_angle(self, i, j, k):
        """bond_angle(i, j, k) --> bond angle defined by atoms i, j, k.
        Angle is calculated using the shortest ji and jk lengths with
        respect to periodic boundary conditions.

        i, j, k  -- atom indices starting at 1

        Return a tuple of (angle, angle_error), both values are in degrees.

        Raises: ValueError if selected atom(s) does not exist
                pdffit.unassignedError when no structure has been loaded
        """
        rv = pdffit2.bond_angle(self._handle, i, j, k)
        return rv


    def blen(self, *args):
        """blen(i, j) --> Show bond length defined by atoms i and j.

        i      -- index of the first atom starting at 1
        j      -- index of the second atom starting at 1

        No return value.  Use bond_length_atoms() to retrieve result.

        Second form:

        blen(a1, a2, lb, ub) --> Show sorted lengths of all a1-a2 bonds.

        a1     -- symbol of the first element in pair or "ALL"
        a2     -- symbol of the second element in pair or "ALL"
        lb     -- lower bond length boundary
        ub     -- upper bond length boundary

        No return value.  Use bond_length_types() to retrieve results.

        Raises: ValueError if selected atom(s) does not exist
                pdffit.unassignedError when no structure has been loaded
        """
        # first form
        if len(args)==2:
            dij, ddij = self.bond_length_atoms(*args[0:2])
            atom_symbols = self.get_atoms()
            ij = (args[0], args[1])
            # indices were already checked in bond_length_atoms call
            assert (0 <= min(ij) - 1) and (max(ij) - 1 < len(atom_symbols))
            symij = ( atom_symbols[ij[0] - 1].upper(),
                      atom_symbols[ij[1] - 1].upper() )
            print(_format_bond_length(dij, ddij, ij, symij), file=output.stdout)
        # second form
        elif len(args)==4:
            a1, a2, lb, ub = args
            try:
                atom_types = self.get_atom_types()
                if isinstance(a1, numbers.Integral):   a1 = atom_types[a1 - 1]
                if isinstance(a2, numbers.Integral):   a2 = atom_types[a2 - 1]
            except IndexError:
                # index of non-existant atom type
                return
            # arguments are OK here, get bond length dictionary
            bld = pdffit2.bond_length_types(self._handle, a1, a2, lb, ub)
            s = "(%s,%s) bond lengths in [%gA,%gA] for current phase :" % \
                    (a1, a2, lb, ub)
            print(s, file=output.stdout)
            atom_symbols = self.get_atoms()
            npts = len(bld['dij'])
            for idx in range(npts):
                dij = bld['dij'][idx]
                ddij = bld['ddij'][idx]
                ij0 = bld['ij0'][idx]
                ij1 = bld['ij1'][idx]
                symij = (atom_symbols[ij0[0]], atom_symbols[ij0[1]])
                s = _format_bond_length(dij, ddij, ij1, symij)
                print(s, file=output.stdout)
            print(file=output.stdout)
            if not bld['dij']:
                print("   *** No pairs found ***", file=output.stdout)
        else:
            emsg = "blen() takes 2 or 4 arguments (%i given)" % len(args)
            raise TypeError(emsg)
        # done
        return


    def bond_length_atoms(self, i, j):
        """bond_length_atoms(i, j) --> shortest distance between atoms i, j.
        Periodic boundary conditions are applied to find the shortest bond.

        i   -- index of the first atom starting at 1
        j   -- index of the second atom starting at 1

        Return a tuple of (distance, distance_error).

        Raises: ValueError if selected atom(s) does not exist
                pdffit.unassignedError when no structure has been loaded.
        """
        rv = pdffit2.bond_length_atoms(self._handle, i, j)
        return rv


    def bond_length_types(self, a1, a2, lb, ub):
        """bond_length_types(a1, a2, lb, ub) --> get all a1-a2 distances.

        a1     -- symbol of the first element in pair or "ALL"
        a2     -- symbol of the second element in pair or "ALL"
        lb     -- lower bond length boundary
        ub     -- upper bond length boundary

        Return a dictionary of distance data containing

        dij  : list of bond lenghts within given bounds
        ddij : list of bond legnth standard deviations
        ij0  : pairs of atom indices starting from 0
        ij1  : pairs of atom indices starting from 1

        Raises: ValueError if selected atom(s) does not exist
                pdffit.unassignedError when no structure has been loaded.
        """
        rv = pdffit2.bond_length_types(self._handle, a1, a2, lb, ub)
        return rv


    def show_scat(self, stype):
        """show_scat(stype) --> Print scattering length for all atoms in
        the current phase.

        stype -- 'X' (xray) or 'N' (neutron).

        Raises: pdffit2.unassignedError if no phase exists
        """
        print(self.get_scat_string(stype), file=output.stdout)
        return


    def get_scat_string(self, stype):
        """get_scat_string(stype) --> Get string with scattering factors
        of all atoms in the current phase.

        stype -- 'X' (xray) or 'N' (neutron).

        Raises:
            pdffit2.unassignedError if no phase exists

        Returns: string with all scattering factors.
        """
        return pdffit2.get_scat_string(self._handle, six.b(stype))


    def get_scat(self, stype, element):
        """get_scat(stype, element) --> Get active scattering factor for
        given element.  If scattering factor has been changed using
        set_scat the result may depend on the active phase.  When no
        phase has been loaded, return the standard value.

        stype   -- 'X' (xray) or 'N' (neutron).
        element -- case-insensitive element symbol such as "Na" or "CL"

        Return float.

        Raises:
            ValueError if element is not known.
        """
        rv = pdffit2.get_scat(self._handle, six.b(stype), element)
        return rv


    def set_scat(self, stype, element, value):
        """set_scat(stype, element, value) --> Set custom scattering factor
        for given element.  The new scattering factor applies only for the
        current phase, in other phases it keeps its default value.

        stype   -- 'X' (xray) or 'N' (neutron).
        element -- case-insensitive element symbol such as "Na" or "CL"
        value   -- new value of scattering factor

        No return value.

        Raises:
            pdffit2.unassignedError if no phase exists.
            ValueError if element is not known.

        See also reset_scat, get_scat.
        """
        pdffit2.set_scat(self._handle, six.b(stype), element, value)
        return


    def reset_scat(self, element):
        """reset_scat(stype, element) --> Reset scattering factors for
        given element to their standard values.  The reset_scat applies
        only for the current phase.

        element -- case-insensitive element symbol such as "Na" or "CL"
        Raises:
            pdffit2.unassignedError if no phase exists
            ValueError if element is not known.
        """
        pdffit2.reset_scat(self._handle, element)
        return


    def num_atoms(self):
        """num_atoms() --> Get number of atoms in current phase.

        Raises: pdffit2.unassignedError if no atoms exist
        """
        return pdffit2.num_atoms(self._handle)


    def num_phases(self):
        """num_phases() --> Number of phases loaded in PdfFit instance.

        Use setphase to bring a specific phase in focus.

        Return integer.
        """
        n = pdffit2.num_phases(self._handle)
        return n


    def num_datasets(self):
        """num_datasets() --> Number of datasets loaded in PdfFit instance.

        Use setdata to bring a specific dataset in focus.

        Return integer.
        """
        n = pdffit2.num_datasets(self._handle)
        return n


    def phase_fractions(self):
        """phase_fractions() --> relative phase fractions for current dataset.
        Convert phase scale factors to relative phase fractions given the
        scattering type of current dataset.

        Return a dictionary of phase fractions with following keys:

        "atom"    -- list of fractions normalized to atom count
        "stdatom" -- errors of atom count fractions
        "cell"    -- list of fractions normalized to unit cell count
        "stdcell" -- errors of unit cell count fractions
        "mass"    -- list of relative weight fractions
        "stdmass" -- errors of relative weight fractions

        Raises: pdffit2.unassignedError if no dataset exists.
        """
        return pdffit2.phase_fractions(self._handle)

    # Begin refinable variables.

    def lat(n):
        """lat(n) --> Get reference to lattice variable n.

        n can be an integer or a string representing the lattice variable.
        1 <==> 'a'
        2 <==> 'b'
        3 <==> 'c'
        4 <==> 'alpha'
        5 <==> 'beta'
        6 <==> 'gamma'
        """
        LatParams = { 'a':1, 'b':2, 'c':3, 'alpha':4, 'beta':5, 'gamma':6 }
        if isinstance(n, six.string_types):
            n = LatParams[n]
        return "lat(%i)" % n
    lat = staticmethod(lat)


    def x(i):
        """x(i) --> Get reference to x-value of atom i."""
        return "x(%i)" % i
    x = staticmethod(x)


    def y(i):
        """y(i) --> Get reference to y-value of atom i."""
        return "y(%i)" % i
    y = staticmethod(y)


    def z(i):
        """z(i) --> Get reference to z-value of atom i."""
        return "z(%i)" % i
    z = staticmethod(z)


    def u11(i):
        """u11(i) --> Get reference to U(1,1) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u11(%i)" % i
    u11 = staticmethod(u11)


    def u22(i):
        """u22(i) --> Get reference to U(2,2) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u22(%i)" % i
    u22 = staticmethod(u22)


    def u33(i):
        """u33(i) --> Get reference to U(3,3) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u33(%i)" % i
    u33 = staticmethod(u33)


    def u12(i):
        """u12(i) --> Get reference to U(1,2) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u12(%i)" % i
    u12 = staticmethod(u12)


    def u13(i):
        """u13(i) --> Get reference to U(1,3) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u13(%i)" % i
    u13 = staticmethod(u13)


    def u23(i):
        """u23(i) --> Get reference to U(2,3) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u23(%i)" % i
    u23 = staticmethod(u23)


    def occ(i):
        """occ(i) --> Get reference to occupancy of atom i."""
        return "occ(%i)" % i
    occ = staticmethod(occ)


    def pscale():
        """pscale() --> Get reference to pscale.

        pscale is the fraction of the total structure that the current phase
        represents.
        """
        return "pscale"
    pscale = staticmethod(pscale)


    def sratio():
        """sratio() --> Get reference to sigma ratio.

        The sigma ratio determines the reduction in the Debye-Waller factor for
        distances below rcut.
        """
        return "sratio"
    sratio = staticmethod(sratio)


    def delta1():
        """delta1() --> Get reference to 1/R peak sharpening factor.
        """
        return "delta1"
    delta1 = staticmethod(delta1)


    def delta2():
        """delta2() --> Reference to (1/R^2) sharpening factor.
        The phenomenological correlation constant in the Debye-Waller factor.
        The (1/R^2) peak sharpening factor.
        """
        return "delta2"
    delta2 = staticmethod(delta2)


    def dscale():
        """dscale() --> Get reference to dscale.

        The data scale factor.
        """
        return "dscale"
    dscale = staticmethod(dscale)


    def qdamp():
        """qdamp() --> Get reference to qdamp.

        Qdamp controls PDF damping due to instrument Q-resolution.
        """
        return "qdamp"
    qdamp = staticmethod(qdamp)


    def qbroad():
        """qbroad() --> Get reference to qbroad.

        Quadratic peak broadening factor.
        """
        return "qbroad"
    qbroad = staticmethod(qbroad)


    def spdiameter():
        """spdiameter() --> Get reference to spdiameter (phase property).

        Diameter value for the spherical particle PDF correction.
        Spherical envelope is not applied when spdiameter equals 0.
        """
        return "spdiameter"
    spdiameter = staticmethod(spdiameter)


    def stepcut():
        """stepcut() --> Get reference to stepcut (phase property).

        stepcut is cutoff radius for empirical step-function PDF envelope.
        stepcut can be used to approximate loss of pair correlations
        in amorphous phase.  stepcut cannot be refined.

        Step cutoff is not applied when stepcut equals 0.
        """
        return "stepcut"
    stepcut = staticmethod(stepcut)


    def rcut():
        """rcut() --> Get reference to rcut.

        rcut is the value of r below which peak sharpening, defined by
        the sigma ratio (sratio), applies.  rcut cannot be refined.
        """
        return "rcut"
    rcut = staticmethod(rcut)


    # End refinable variables.

    def __init__(self):

        self.stru_files = []
        self.data_files = []

        self._handle = pdffit2.create()
        self.intro()
        return


    def __getRef(self, var_string):
        """Return the actual reference to the variable in the var_string.

        This function must be called before trying to actually reference an
        internal variable. See the constrain method for an example.

        Raises:
            pdffit2.unassignedError if variable is not yet assigned
            ValueError if variable index does not exist (e.g. lat(7))
        """
        var_string = _convertCallable(var_string)
        arg_int = None
        try:
            method_string, arg_string = var_string.split("(")
            method_string = method_string.strip()
            arg_int = int(arg_string.strip(")").strip())
        except ValueError: #There is no arg_string
            method_string = var_string.strip()

        f = getattr(pdffit2, method_string)
        if arg_int is None:
            retval = f(self._handle)
        else:
            retval = f(self._handle, arg_int)
        return retval


# End of class PdfFit


# End of file
