from GSASII import GSASIIElem

def test_get_xsection():
    xsection = GSASIIElem.GetXsectionCoeff('Fe')
    assert len(xsection) > 0
