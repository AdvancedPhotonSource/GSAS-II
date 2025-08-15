"""
test_nistlat.py
===============
Tests the two NIST*LATTICE binaries, LATTIC and convcell using
their Python wrappers in module nistlat.py
"""

import importlib.util
import os
import sys
import tempfile

import numpy.testing as npt

G2loc = None
try:
    G2loc = importlib.util.find_spec("GSASII.nistlat")
except ModuleNotFoundError:
    print("ModuleNotFound for GSASII.nistlat")

home = os.path.dirname(__file__)
if G2loc is None:  # fixup path if GSASII not installed into Python
    print("GSAS-II not installed in Python; Hacking sys.path")
    sys.path.append(os.path.dirname(home))

import GSASII.GSASIIlattice as G2lat
from GSASII import nistlat

V = lambda cell: float(G2lat.calc_V(G2lat.cell2A(cell)))

work = tempfile.gettempdir()


def test_CellSymSearch():
    "test_CellSymSearch(): nistlat.CellSymSearch"
    print("test_CellSymSearch(): nistlat.CellSymSearch")
    cell = [14.259, 22.539, 8.741, 90.0, 114.1, 90.0]
    center = "C"
    # test CellSymSearch with out sub- & supercells (default mode)
    msg = "CellSymSearch #1 reduced cell"
    res = nistlat.CellSymSearch(cell, center)
    npt.assert_allclose(
        res[0][1][0],
        [8.741, 13.3353, 13.3393, 115.369, 102.638, 102.61],
        rtol=0.001,
        err_msg=msg,
    )
    assert res[0][1][1] == "P", msg

    msg = "CellSymSearch #1 centered cell"
    npt.assert_allclose(
        res[0][2][0],
        [22.5417, 22.5417, 8.741, 90.0, 90.0, 120.0],
        rtol=0.001,
        err_msg=msg,
    )
    assert res[0][2][1] == "R" and res[0][2][2] == "H", msg

    # now test CellSymSearch in mode=3 (w/sub- & supercells)
    msg = "CellSymSearch #2 centered cell w/sub & supercells"
    res = nistlat.CellSymSearch(cell, center, mode=3)
    vPrim = [V(r[1][0]) / V(cell) for r in res]  # cell Vol of primitive cell
    npt.assert_allclose(
        vPrim,
        [
            0.5,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
        ],
        rtol=0.001,
        err_msg=msg,
    )

    vCent = [V(r[2][0]) / V(cell) for r in res]  # cell Vol of full cell
    npt.assert_allclose(
        vCent,
        [1.5, 1.0, 1.0, 4.0, 3.0, 1.0, 4.0, 4.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75],
        rtol=0.001,
        err_msg=msg,
    )
    print("OK")


def test_CompareCell():
    "test_CompareCell(): nistlat.CompareCell"
    print("test_CompareCell(): nistlat.CompareCell")
    cell1 = (5.03461, 5.03461, 13.74753, 90, 90, 120)
    center1 = "R"
    cell2 = (7.40242, 5.03461, 5.42665, 90, 84.14, 90)
    center2 = "P"
    center2 = "I"
    tolerance = 3 * [0.1] + 3 * [0.5]
    res = nistlat.CompareCell(cell1, center1, cell2, center2, tolerance)
    msg = "CompareCell number of matches"
    assert len(res) == 6, msg
    for i, r in enumerate(res):
        msg = f"CompareCell result #{i}"
        # check the forward and reverse cell transforms
        npt.assert_allclose(
            cell1, G2lat.TransformCell(cell2, r[5])[:6], rtol=0.001, err_msg=msg
        )
        npt.assert_allclose(
            cell2, G2lat.TransformCell(cell1, r[4])[:6], rtol=0.001, err_msg=msg
        )
    print("OK")


def test_ConvCell():
    "test_ConvCell(): nistlat.ConvCell"
    print("test_ConvCell(): nistlat.ConvCell")
    # convert a rhombohedral cell to hexagonal setting
    cellin = [
        5.0,
        5.0,
        5.0,
        85.0,
        85.0,
        85.0,
    ]
    cellout = nistlat.ConvCell(cellin)
    msg = "ConvCell rhomb->hex cell"
    assert cellout[1] == "R" and cellout[2] == "H", msg
    npt.assert_allclose(
        cellout[0], [6.7559, 6.7559, 9.3847, 90.0, 90.0, 120.0], rtol=0.001, err_msg=msg
    )
    # a hexagonal cell is unchanged
    msg = "ConvCell unchanged hex cell"
    cellin = cellout[0]
    cellout = nistlat.ConvCell(cellin)
    assert cellout[1] == "P", msg
    npt.assert_allclose(cellout[0], cellin, rtol=0.001, err_msg=msg)
    print("OK")


def test_ReduceCell():
    "test_ReduceCell(): nistlat.ReduceCell"
    print("test_ReduceCell(): nistlat.ReduceCell")
    res = nistlat.ReduceCell("I", [3, 3, 5, 90, 90, 90])  # body-center tetragonal
    rcell = res["output"][0][1]
    msg = "ReduceCell body-center tetragonal cell"
    npt.assert_allclose(
        rcell, [3.0, 3.0, 3.2787, 117.226, 117.226, 90.0], rtol=0.001, err_msg=msg
    )
    ocell = nistlat.ConvCell(rcell)
    assert ocell[1] == "I", msg
    npt.assert_allclose(
        ocell[0], [3.0, 3.0, 5.0, 90.0, 90.0, 90.0], rtol=0.001, err_msg=msg
    )
    print("OK")


if __name__ == "__main__":
    test_ConvCell()
    test_ReduceCell()
    test_CellSymSearch()
    test_CompareCell()
