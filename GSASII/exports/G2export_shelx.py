"""Classes in :mod:`~GSASII.exports.G2export_shelx` follow:"""

import os.path

from .. import GSASIIfiles as G2fil
from .. import GSASIIspc as G2spc


class ExportPhaseShelx(G2fil.ExportBaseclass):
    """Used to create a SHELX .ins file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    """

    def __init__(self, G2frame):
        super(self.__class__, self).__init__(  # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName="SHELX .ins",
            extension=".ins",
            longFormatName="Export phase as SHELX .ins file",
        )
        self.exporttype = ["phase"]
        self.multiple = True

    def Exporter(self, event=None):
        """Export as a SHELX .ins file"""
        import re

        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()
        if self.ExportSelect():
            return  # set export parameters;
        filename = self.filename
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam]  # pointer to current phase info
            i = self.Phases[phasenam]["pId"]
            if (
                len(self.phasenam) > 1
            ):  # if more than one filename is included, add a phase #
                self.filename = (
                    os.path.splitext(filename)[1] + "_" + str(i) + self.extension
                )
            self.OpenFile(filename)
            # title line
            self.Write(
                "TITL from "
                + str(self.G2frame.GSASprojectfile)
                + ", phase "
                + str(phasenam)
            )
            # get & write cell parameters
            cell, sig = self.GetCell(phasenam)
            self.Write(
                "CELL 0.5 {:.5f} {:.5f} {:.5f} {:.3f} {:.3f} {:.3f}".format(*cell[:6])
            )
            # Shelx lattice number
            lattnum = {"P": 1, "I": 2, "R": 2, "F": 3, "A": 4, "B": 5, "C": 6}.get(
                phasedict["General"]["SGData"]["SGLatt"], 0
            )
            if not phasedict["General"]["SGData"]["SGInv"]:
                lattnum *= -1
            self.Write("LATT " + str(lattnum))
            # generate symmetry operations not including centering and center of symmetry
            for Opr in phasedict["General"]["SGData"]["SGOps"]:
                sym = G2spc.MT2text(Opr).lower().replace(" ,", ", ")
                self.Write("SYMM " + G2fil.trim(sym))
            # scan through atom types, count the number of times that each element occurs
            AtomsList = self.GetAtoms(phasenam)
            maxmult = 0
            elemtypes = {}
            for lbl, typ, mult, xyz, td in AtomsList:
                maxmult = max(maxmult, mult)
                typ = re.search("[A-Za-z]+", typ).group(0)
                typ = typ[0:2]
                typ = typ[0].upper() + typ[1:].lower()
                if elemtypes.get(typ) is None:
                    elemtypes[typ] = 1
                else:
                    elemtypes[typ] += 1
            # create scattering factor record
            s = "SFAC"
            elemlist = sorted(elemtypes)
            for elem in elemlist:
                s += " " + elem
            self.Write(s)
            # handle atom records
            count = {}
            for lbl, typ, mult, xyz, td in AtomsList:
                typ = re.search("[A-Za-z]+", typ).group(0)
                typ = typ[0:2]
                typ = typ[0].upper() + typ[1:].lower()
                if count.get(typ) is None:
                    count[typ] = 0
                else:
                    count[typ] += 1
                # make a unique <=4 character label, if possible
                if elemtypes[typ] <= 99:
                    lbl = f"{typ:s}{count[typ]:d}"
                else:  # more than 99 atoms, use hexadecimal notation
                    lbl = typ + f"{count[typ]:X}"[-2:]
                sfnum = (
                    elemlist.index(typ) + 1
                )  # element number in scattering factor list
                l = lbl + " " + str(sfnum)
                l += " {:.5f} {:.5f} {:.5f}".format(*[x[0] for x in xyz[:3]])
                if mult == maxmult:
                    m = 1
                else:
                    m = 1.0 * mult / maxmult
                if xyz[3][0] == 1:  # frac
                    occ = 10 + m
                else:
                    occ = m * xyz[3][0]
                l += f" {occ:.3f}"
                for val, _sig in td:
                    l += f" {val:.3f}"
                self.Write(l)
            self.Write("END")
            self.CloseFile()
            print("Phase " + phasenam + " written to file " + self.fullpath)
