"""
test_diffax.py
==============
Performs a DIFFaX computation using the GSAS-II interface to the
pydiffax.so [.pyd] module
"""

import os
import sys
import tempfile

import numpy as np
import numpy.testing as npt

home = os.path.dirname(__file__)
work = tempfile.gettempdir()

import importlib.util

G2loc = None
try:
    G2loc = importlib.util.find_spec("GSASII.GSASIIpwd")
except ModuleNotFoundError:
    print("ModuleNotFound for GSASII.GSASIIpwd")

if G2loc is None:  # fixup path if GSASII not installed into Python
    print("GSAS-II not installed in Python; Hacking sys.path")
    sys.path.append(os.path.dirname(home))

# import GSASII.GSASIIscriptable as G2sc  # sets up access to binaries
import GSASII.GSASIIpwd as G2pwd


def test_diffax():
    "tests DIFFaX"
    print("test_diffax(): test a small DIFFaX computation")
    # define input needed for computation
    layerinfo = {
        "Layers": [
            {
                "SameAs": "",
                "Name": "layer 1",
                "Symm": "-1",
                "Atoms": [["C1", "C", -1.0 / 3.0, -3.0 / 18.0, -1 / 8, 1.0, 0.01]],
            },
            {
                "SameAs": "",
                "Name": "layer 2",
                "Symm": "-1",
                "Atoms": [["C1", "C", 1.0 / 3.0, 3.0 / 18.0, -1 / 8, 1.0, 0.01]],
            },
        ],
        "AtInfo": {
            "C": {
                "Isotopes": {
                    "13": {"SA": 0.00137, "Mass": 13.003, "SL": [0.619, 0]},
                    "12": {"SA": 0.00353, "Mass": 12.0, "SL": [0.6654, 0]},
                    "Nat. Abund.": {"SA": 0.0035, "Mass": 12.011, "SL": [0.6648, 0]},
                },
                "Vdrad": 1.95,
                "Color": (144, 144, 144),
                "Symbol": "C",
                "Lande g": 2.0,
                "Drad": 1.12,
                "Mass": 12.011,
                "Arad": 0.92,
                "Z": 6,
                "Hbrad": 0,
            }
        },
        "allowedTrans": [["1", "1"], ["1", "2"], ["2", "1"], ["2", "2"]],
        "seqCodes": ("TransP;0;0", [0.0, 1.0], 10),
        "SymTrans": True,
        "Stacking": ["recursive", "infinite", ""],
        "Laue": "6/mmm",
        "Cell": [False, 2.522, 2.522, 2.059, 90.0, 90.0, 120.0, 11.341673551466423],
        "Width": [[1.0, 1.0], [False, False]],
        "Sadp": {},
        "seqResults": [],
        "Toler": 0.01,
        "Transitions": [
            [
                [0.7, 2.0 / 3.0, 1.0 / 3.0, 1.0, "", False],
                [0.3, 0.0, 0.0, 1.0, "", False],
            ],
            [
                [0.3, 0.0, 0.0, 1.0, "", False],
                [0.7, -2.0 / 3.0, -1.0 / 3.0, 1.0, "", False],
            ],
        ],
        "selInst": "Gaussian",
    }
    bkg = [
        ["chebyschev", True, 3, 50.0, 0.0, 0.0],
        {
            "peaksList": [],
            "debyeTerms": [],
            "nPeaks": 0,
            "background PWDR": ["", -1.0],
            "nDebye": 0,
        },
    ]
    profile = [
        np.linspace(10, 150, 7001),
        np.zeros(7001),
        np.zeros(7001),
        np.zeros(7001),
        np.zeros(7001),
        np.zeros(7001),
    ]
    inst = {
        "I(L2)/I(L1)": [0.5, 0.5, False],
        "SH/L": [0.002, 0.002, False],
        "V": [-2.0, -2.0, False],
        "Lam2": [1.5443, 1.5443, False],
        "Source": ["CuKa", "?"],
        "Zero": [0.0, 0.0, False],
        "Lam1": [1.5405, 1.5405, False],
        "U": [2.0, 2.0, False],
        "W": [5.0, 5.0, False],
        "Azimuth": [0.0, 0.0, False],
        "Y": [0.0, 0.0, False],
        "X": [0.0, 0.0, False],
        "Type": ["PXC", "PXC", False],
        "Bank": [1.0, 1.0, False],
        "Polariz.": [0.7, 0.7, False],
    }

    G2pwd.CalcStackingPWDR(layerinfo, 1000.0, bkg, [10.0, 150.0], inst, profile, False)
    # import matplotlib.pyplot as plt
    # plt.plot(profile[0],profile[3])
    # plt.show()
    assert abs(profile[3][:7000].max() - 5230.06) < 0.1, "Max value off"
    assert abs(profile[3][:7000].min() - 50.03) < 0.1, "Min value off"
    msg = "deviations in 42-43 deg section"
    npt.assert_allclose(
        profile[3][1600:1651],
        np.array(
            [
                200.48310247,
                200.63761559,
                200.8637708,
                201.15861368,
                201.51951493,
                201.94413305,
                202.43038195,
                202.97640301,
                203.58054064,
                204.24132092,
                204.95743283,
                205.72771164,
                206.55112422,
                207.42675591,
                208.35379877,
                209.33154098,
                210.35935728,
                211.43670024,
                212.56309224,
                213.73811812,
                214.96141832,
                216.23268245,
                217.55164331,
                218.91807113,
                220.33176811,
                221.79256313,
                223.30030665,
                224.85486568,
                226.45611886,
                228.10395155,
                229.79825094,
                231.53890116,
                233.32577831,
                235.15874548,
                237.03764759,
                238.96230622,
                240.93251424,
                242.94803037,
                245.00857357,
                247.11381723,
                249.26338327,
                251.45683604,
                253.69367607,
                255.97333369,
                258.29516247,
                260.6584326,
                263.06232408,
                265.50591998,
                267.98819954,
                270.5080314,
                273.06416684,
            ]
        ),
        rtol=0.0001,
        err_msg=msg,
    )


if __name__ == "__main__":
    # import time
    # start = time.time()
    test_diffax()
    # print('elapsed=',time.time()-start)
