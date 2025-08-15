import pickle

from .. import GSASIIobj as G2obj


class G2_ReaderClass(G2obj.ImportImage):
    """Routine to read an image that has been pickled in Python. Images
    in this format are created by the "Sum image data" command. At least for
    now, only one image is permitted per file.
    """

    def __init__(self):
        super(self.__class__, self).__init__(  # fancy way to self-reference
            extensionlist=(".G2img",),
            strictExtension=True,
            formatName="GSAS-II image",
            longFormatName="pickled image from GSAS-II",
        )

    def ContentsValidator(self, filename):
        """test by trying to unpickle (should be quick)"""
        try:
            fp = open(filename, "rb")  # noqa: SIM115
            pickle.load(fp)
            fp.close()
        except:  # noqa: E722
            return False
        return True

    def Reader(self, filename, ParentFrame=None, **unused):  # noqa: ARG002
        """Read using pickle"""
        Fp = open(filename, "rb")  # noqa: SIM115
        self.Comments, self.Data, self.Npix, self.Image = pickle.load(Fp)
        Fp.close()
        self.LoadImage(ParentFrame, filename)
        return True
