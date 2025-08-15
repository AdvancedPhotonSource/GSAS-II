# output from GSASIIspc computed on platform darwin on 2013-06-02
import numpy as np

array = np.array
float32 = np.float32
# testing 255 space groups (25 dups/non-standard)
SGdat = {
    "p 4/n b m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/n b m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4 c 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 c 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "r -3 m": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "",
        "SGLatt": "R",
        "SpGrp": "R -3 m",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 42 n m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 42 n m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "a b a 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "A",
        "SpGrp": "A b a 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 42/m b c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/m b c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p m n 21": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P m n 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "i 4/m c m ": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 4/m c m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42/m c m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/m c m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p m -3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m -3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p b a 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P b a 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i b a m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I b a m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 21/m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 21/m",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 41": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 41",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.75], dtype=float32),
            ],
        ],
    },
    "p 42": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 42",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 43": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 43",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.25], dtype=float32),
            ],
        ],
    },
    "f 4 3 2": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F 4 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4 21 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 21 m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 63/m c m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 63/m c m",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 2 3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 2 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i a 3": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I a 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 3 2 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 3 2 1",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i a -3 d": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I a -3 d",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
        ],
    },
    "p a -3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P a -3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "c 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "y",
        "SGLatt": "C",
        "SpGrp": "C 2",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "r 3 2 h": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "",
        "SGLatt": "R",
        "SpGrp": "R 3 2 h",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "x z",
        "SGLatt": "P",
        "SpGrp": "P c",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "f 2 2 2": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F 2 2 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 6 c c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 6 c c",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "r 3 2 r": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "R 3 2 r",
        "SGLaue": "3mR",
        "SGSys": "rhombohedral",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 62 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 62 2 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
        ],
    },
    "i 41 c d": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I 41 c d",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.25], dtype=float32),
            ],
        ],
    },
    "f m m 2": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "z",
        "SGLatt": "F",
        "SpGrp": "F m m 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p m m 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P m m 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "x z",
        "SGLatt": "P",
        "SpGrp": "P m",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 4 2 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 4 2 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 31 2 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 31 2 1",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
        ],
    },
    "i -4": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I -4",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "y",
        "SGLatt": "P",
        "SpGrp": "P 2",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 21 21 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 21 21 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c 1 2/c 1": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C 1 2/c 1",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "i b a 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I b a 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p b a m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P b a m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p b a n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P b a n",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p m -3 n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m -3 n",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "i b c a": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I b c a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42 21 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42 21 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "f m m m": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F m m m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 41": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I 41",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
        ],
    },
    "p 6": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 6",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 3",
        "SGLaue": "3",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p m m n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m m n",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p m m m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m m m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "r 3 m h": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "z",
        "SGLatt": "R",
        "SpGrp": "R 3 m h",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c m c 21": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "z",
        "SGLatt": "C",
        "SpGrp": "C m c 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "xyz",
        "SGLatt": "P",
        "SpGrp": "P 1",
        "SGLaue": "-1",
        "SGSys": "triclinic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ]
        ],
    },
    "i 4": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I 4",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 4": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 4",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 42 b c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 42 b c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p m m a": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m m a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i -4 m 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I -4 m 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4 21 c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 21 c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4 2 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 4/m c c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/m c c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p -6 2 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -6 2 m",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 1 2/m 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 1 2/m 1",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -6 2 c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -6 2 c",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 6 m m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 6 m m",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c c": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "x z",
        "SGLatt": "C",
        "SpGrp": "C c",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 43 3 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 43 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "a 2 2 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "A",
        "SpGrp": "A 2 2 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -3",
        "SGLaue": "3",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -1",
        "SGLaue": "-1",
        "SGSys": "triclinic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ]
        ],
    },
    "f d d 2": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "z",
        "SGLatt": "F",
        "SpGrp": "F d d 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 62": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 62",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
        ],
    },
    "c m m a": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C m m a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p -3 c 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -3 c 1",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "c m c m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C m c m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "c m m m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C m m m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c m c a": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C m c a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "i a -3": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I a -3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "i m a 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I m a 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 63/m m c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 63/m m c",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4 3 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 41/a c d": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 41/a c d",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
        ],
    },
    "p -4 2 c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 2 c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4 n c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 4 n c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4/m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/m",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 4/n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/n",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 21/c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 21/c",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "r -3 c": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "",
        "SGLatt": "R",
        "SpGrp": "R -3 c",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4/n m m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/n m m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 4/m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 4/m",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 3 m 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 3 m 1",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 63/m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 63/m",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 6 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 6 2 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 2/m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 2/m",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "f d 3": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F d 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.75, 0.25], dtype=float32),
            ],
        ],
    },
    "i 41/a m d": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 41/a m d",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
        ],
    },
    "p 4/n c c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/n c c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "i m m a": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I m m a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 4 b m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 4 b m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 2/c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 2/c",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p -6 m 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -6 m 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p n n 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P n n 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 31 1 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 31 1 2",
        "SGLaue": "31m",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
        ],
    },
    "f -4 3 c": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F -4 3 c",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i m -3": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I m -3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "f -4 3 m": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F -4 3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 21 3": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 21 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42/m m c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/m m c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 65 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 65 2 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.83333331], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.16666667], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.83333331], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.16666667], dtype=float32),
            ],
        ],
    },
    "p 4/m n c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/m n c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "c 2/m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C 2/m",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "f d d d": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F d d d",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.5], dtype=float32),
            ],
        ],
    },
    "c m m 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "z",
        "SGLatt": "C",
        "SpGrp": "C m m 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 43 21 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 43 21 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -3 1 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -3 1 m",
        "SGLaue": "31m",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 2 2 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 2 2 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 42/n b c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/n b c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "i 4 3 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 4 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 41 3 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 41 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42/n m c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/n m c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 64 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 64 2 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
        ],
    },
    "p c a 21": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P c a 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "f d -3 c": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F d -3 c",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.75, 0.75], dtype=float32),
            ],
        ],
    },
    "p n a 21": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P n a 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p -4 n 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 n 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42/n n m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/n n m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "f d -3 m": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F d -3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.75, 0.25], dtype=float32),
            ],
        ],
    },
    "r -3 c r": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "R -3 c r",
        "SGLaue": "3mR",
        "SGSys": "rhombohedral",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 63 m c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 63 m c",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4/m b m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/m b m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 2 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 2 2 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 63 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 63 2 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 6/m m m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 6/m m m",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p c c n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P c c n",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p c c m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P c c m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p m n a": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m n a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "f 41 3 2": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F 41 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.25, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "r -3 r": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "R -3 r",
        "SGLaue": "3R",
        "SGSys": "rhombohedral",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 1 1 2/m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 1 1 2/m",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "c",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 64": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 64",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
        ],
    },
    "p c c a": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P c c a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "f m -3": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F m -3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -6": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -6",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i m m m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I m m m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4 2 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 2 m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 21 3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 21 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4 m m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 4 m m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4 m 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 m 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c 2/c": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C 2/c",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42 3 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 6/m c c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 6/m c c",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "f m 3": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F m 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p n n a": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P n n a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i -4 3 d": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I -4 3 d",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p n n n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P n n n",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p n n m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P n n m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i -4 3 m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I -4 3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 65": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 65",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.83333331], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.16666667], dtype=float32),
            ],
        ],
    },
    "r 3 r": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "R 3 r",
        "SGLaue": "3R",
        "SGSys": "rhombohedral",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 2/m 1 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 2/m 1 1",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": True,
        "SGUniq": "a",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 41/a": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 41/a",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
        ],
    },
    "p 63 c m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 63 c m",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c 1 2 1": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "y",
        "SGLatt": "C",
        "SpGrp": "C 1 2 1",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p b c n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P b c n",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p b c m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P b c m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "a m m 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "A",
        "SpGrp": "A m m 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i m -3 m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I m -3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 4 m m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I 4 m m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 61 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 61 2 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.16666667], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.83333331], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.16666667], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.83333331], dtype=float32),
            ],
        ],
    },
    "i m m 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I m m 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 42/n c m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/n c m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p b c a": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P b c a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4 21 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4 21 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 4/n n c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/n n c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "f m -3 m": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F m -3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 4/m m m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 4/m m m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "f m -3 c": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F m -3 c",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p n -3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P n -3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p c c 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P c c 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 41 3 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 41 3 2",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
                    dtype=float32,
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]],
                    dtype=float32,
                ),
                array([0.75, 0.75, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42 m c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 42 m c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 4 c c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 4 c c",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p m -3 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m -3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 32 1 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 32 1 2",
        "SGLaue": "31m",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
        ],
    },
    "p 32 1 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 32 1 1",
        "SGLaue": "3",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
        ],
    },
    "r -3 m r": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "R -3 m r",
        "SGLaue": "3mR",
        "SGSys": "rhombohedral",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 3 c 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 3 c 1",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 2 2 21": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 2 2 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 63": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 63",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p m 3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P m 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 42/m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/m",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p m c 21": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P m c 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42/n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/n",
        "SGLaue": "4/m",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "a m a 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "A",
        "SpGrp": "A m a 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 6/m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 6/m",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -6 c 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -6 c 2",
        "SGLaue": "6/mmm",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i -4 c 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I -4 c 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "F -1": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F -1",
        "SGLaue": "-1",
        "SGSys": "triclinic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ]
        ],
    },
    "p 3 1 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 3 1 m",
        "SGLaue": "31m",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c c c 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "z",
        "SGLatt": "C",
        "SpGrp": "C c c 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i m 3": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I m 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4 3 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -4 3 n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 3 n",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p -3 1 c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -3 1 c",
        "SGLaue": "31m",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "r 3 m": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "z",
        "SGLatt": "R",
        "SpGrp": "R 3 m",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 21": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "y",
        "SGLatt": "P",
        "SpGrp": "P 21",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "r -3": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "",
        "SGLatt": "R",
        "SpGrp": "R -3",
        "SGLaue": "3",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "x z",
        "SGLatt": "C",
        "SpGrp": "C m",
        "SGLaue": "2/m",
        "SGSys": "monoclinic",
        "SGInv": False,
        "SGUniq": "b",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 32 2 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 32 2 1",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
        ],
    },
    "i 21 21 21": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 21 21 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42 2 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "i -4 2 m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I -4 2 m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 65 1 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 65 1 1",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.83333331], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.16666667], dtype=float32),
            ],
        ],
    },
    "p 61": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 61",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.16666667], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.83333331], dtype=float32),
            ],
        ],
    },
    "i 2 3": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 2 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i -4 2 d": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I -4 2 d",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
        ],
    },
    "p a 3": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P a 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "f 2 3": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F 2 3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 4 c m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I 4 c m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "r 3 c": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "z",
        "SGLatt": "R",
        "SpGrp": "R 3 c",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p n m a": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P n m a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "r 3 c r": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "R 3 c r",
        "SGLaue": "3mR",
        "SGSys": "rhombohedral",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p n c 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P n c 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c 2 2 21": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C 2 2 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "r 3 m r": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "R 3 m r",
        "SGLaue": "3mR",
        "SGSys": "rhombohedral",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 43 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 43 2 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.25], dtype=float32),
            ],
        ],
    },
    "r 3 2": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "",
        "SGLatt": "R",
        "SpGrp": "R 3 2",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p m a 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P m a 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 4/m m m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 4/m m m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c c c a": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C c c a",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 41 m d": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "I",
        "SpGrp": "I 41 m d",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
        ],
    },
    "c c c m": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C c c m",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 41 21 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 41 21 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 31": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 31",
        "SGLaue": "3",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
        ],
    },
    "p 32": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 32",
        "SGLaue": "3",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.66666669], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.33333334], dtype=float32),
            ],
        ],
    },
    "p 42/m n m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 42/m n m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 3 1 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 3 1 2",
        "SGLaue": "31m",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "i 41 2 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        "SGPolax": "",
        "SGLatt": "I",
        "SpGrp": "I 41 2 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p -3 m 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -3 m 1",
        "SGLaue": "3m1",
        "SGSys": "trigonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "a b m 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5]]),
        "SGPolax": "z",
        "SGLatt": "A",
        "SpGrp": "A b m 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p n -3 n": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P n -3 n",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "r 3": {
        "SGCen": array(
            [
                [0.0, 0.0, 0.0],
                [0.33333333, 0.66666667, 0.66666667],
                [0.66666667, 0.33333333, 0.33333333],
            ]
        ),
        "SGPolax": "z",
        "SGLatt": "R",
        "SpGrp": "R 3",
        "SGLaue": "3",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "c 2 2 2": {
        "SGCen": array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        "SGPolax": "",
        "SGLatt": "C",
        "SpGrp": "C 2 2 2",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p n -3 m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P n -3 m",
        "SGLaue": "m3m",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
        ],
    },
    "p 42 c m": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 42 c m",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 6/m 1 1": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 6/m 1 1",
        "SGLaue": "6/m",
        "SGSys": "hexagonal",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
        ],
    },
    "p 21 21 21": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 21 21 21",
        "SGLaue": "mmm",
        "SGSys": "orthorhombic",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.5, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "f d -3": {
        "SGCen": array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        ),
        "SGPolax": "",
        "SGLatt": "F",
        "SpGrp": "F d -3",
        "SGLaue": "m3",
        "SGSys": "cubic",
        "SGInv": True,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.25, 0.25, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.0, 0.25, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.5, 0.75], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.75, 0.25, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], dtype=float32
                ),
                array([0.5, 0.25, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.75, 0.5, 0.25], dtype=float32),
            ],
            [
                array(
                    [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float32
                ),
                array([0.25, 0.75, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.75, 0.25], dtype=float32),
            ],
        ],
    },
    "p -4 b 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P -4 b 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.5, 0.5, 0.0], dtype=float32),
            ],
        ],
    },
    "p 3 1 c": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "z",
        "SGLatt": "P",
        "SpGrp": "P 3 1 c",
        "SGLaue": "31m",
        "SGSys": "trigonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
        ],
    },
    "p 41 2 2": {
        "SGCen": array([[0, 0, 0]]),
        "SGPolax": "",
        "SGLatt": "P",
        "SpGrp": "P 41 2 2",
        "SGLaue": "4/mmm",
        "SGSys": "tetragonal",
        "SGInv": False,
        "SGUniq": "",
        "SGOps": [
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.75], dtype=float32),
            ],
            [
                array(
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.0], dtype=float32),
            ],
            [
                array(
                    [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                    dtype=float32,
                ),
                array([0.0, 0.0, 0.25], dtype=float32),
            ],
            [
                array(
                    [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.5], dtype=float32),
            ],
            [
                array(
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], dtype=float32
                ),
                array([0.0, 0.0, 0.75], dtype=float32),
            ],
        ],
    },
}
SGlist = {
    "p 4/n b m": [
        " Space Group: P 4/n b m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,    Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,    Z \t",
        " ( 5)    -X ,1/2+Y ,    Z \t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7) 1/2+X ,   -Y ,    Z \t( 8) 1/2+Y ,1/2+X ,    Z \t",
        " ",
    ],
    "p -4 c 2": [
        " Space Group: P -4 c 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)     Y ,    X ,1/2-Z \t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)    -Y ,   -X ,1/2-Z \t",
        " ",
    ],
    "r -3 m": [
        " Space Group: R -3 m",
        " The lattice is centrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 36",
        " The Laue symmetry is 3m1",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,    Z \t\t",
        " ( 5)    -Y ,   -X ,    Z \t\t( 6)     X ,   X-Y,    Z \t\t",
        " ",
    ],
    "p 42 n m": [
        " Space Group: P 42 n m",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,1/2+X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4) 1/2+Y ,1/2-X ,1/2+Z \t",
        " ( 5) 1/2-X ,1/2+Y ,1/2+Z \t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7) 1/2+X ,1/2-Y ,1/2+Z \t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "a b a 2": [
        " Space Group: A b a 2",
        " The lattice is noncentrosymmetric A-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3) 1/2+X ,1/2-Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 42/m b c": [
        " Space Group: P 42/m b c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ( 5) 1/2-X ,1/2+Y ,    Z \t( 6) 1/2-Y ,1/2-X ,1/2+Z \t",
        " ( 7) 1/2+X ,1/2-Y ,    Z \t( 8) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ",
    ],
    "p m n 21": [
        " Space Group: P m n 21",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3) 1/2+X ,   -Y ,1/2+Z \t( 4) 1/2-X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "i 4/m c m ": [
        " Space Group: I 4/m c m",
        " The lattice is centrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 32",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p 42/m c m": [
        " Space Group: P 42/m c m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p m -3": [
        " Space Group: P m -3",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " ",
    ],
    "p b a 2": [
        " Space Group: P b a 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3) 1/2+X ,1/2-Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "i b a m": [
        " Space Group: I b a m",
        " The lattice is centrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3) 1/2+X ,1/2-Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 21/m": [
        " Space Group: P 21/m",
        " The lattice is centrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,   -Z \t",
        " ",
    ],
    "p 41": [
        " Space Group: P 41",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 4/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/4+Z \t",
        " ( 3)    -X ,   -Y ,1/2+Z \t( 4)     Y ,   -X ,3/4+Z \t",
        " ",
    ],
    "p 42": [
        " Space Group: P 42",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 4/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ",
    ],
    "p 43": [
        " Space Group: P 43",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 4/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,3/4+Z \t",
        " ( 3)    -X ,   -Y ,1/2+Z \t( 4)     Y ,   -X ,1/4+Z \t",
        " ",
    ],
    "f 4 3 2": [
        " Space Group: F 4 3 2",
        " The lattice is noncentrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)    -Y ,    X ,    Z \t\t",
        " ( 5)     Z ,   -Y ,    X \t\t( 6)     X ,    Z ,   -Y \t\t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13)     Y ,   -X ,    Z \t\t(14)     Z ,    Y ,   -X \t\t",
        " (15)    -X ,    Z ,    Y \t\t(16)    -X ,   -Z ,   -Y \t\t",
        " (17)    -Y ,   -X ,   -Z \t\t(18)    -Z ,   -Y ,   -X \t\t",
        " (19)     Y ,   -Z ,   -X \t\t(20)     Y ,    X ,   -Z \t\t",
        " (21)    -Z ,    Y ,    X \t\t(22)     X ,   -Z ,    Y \t\t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "p -4 21 m": [
        " Space Group: P -4 21 m",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,   -Z \t( 6) 1/2+Y ,1/2+X ,    Z \t",
        " ( 7) 1/2+X ,1/2-Y ,   -Z \t( 8) 1/2-Y ,1/2-X ,    Z \t",
        " ",
    ],
    "p 63/m c m": [
        " Space Group: P 63/m c m",
        " The lattice is centrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is 6/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/2+Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,1/2+Z \t",
        " ( 7)    Y-X,    Y ,1/2+Z \t( 8)    -X ,   Y-X,    Z \t\t",
        " ( 9)    -Y ,   -X ,1/2+Z \t(10)    X-Y,   -Y ,    Z \t\t",
        " (11)     X ,   X-Y,1/2+Z \t(12)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p 2 3": [
        " Space Group: P 2 3",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is m3",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,   -Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,   -Y \t\t( 6)    -Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,   -X ,    Y \t\t( 8)     Y ,   -Z ,   -X \t\t",
        " ( 9)    -Y ,    Z ,   -X \t\t(10)    -X ,   -Y ,    Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "i a 3": [
        " Space Group: I a 3",
        " The lattice is centrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,    Y ,1/2-Z \t",
        " ( 5) 1/2-Z ,1/2+X ,    Y \t( 6)     Y ,1/2-Z ,1/2+X \t",
        " ( 7)    -Z ,1/2+X ,1/2-Y \t( 8) 1/2-Y ,   -Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/2-Z ,   -X \t(10)    -X ,1/2+Y ,1/2-Z \t",
        " (11) 1/2-Z ,   -X ,1/2+Y \t(12) 1/2+X ,1/2-Y ,   -Z \t",
        " ",
    ],
    "p 3 2 1": [
        " Space Group: P 3 2 1",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3m1",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,   Y-X,   -Z \t\t( 6)    X-Y,   -Y ,   -Z \t\t",
        " ",
    ],
    "i a -3 d": [
        " Space Group: I a -3 d",
        " The lattice is centrosymmetric I-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,    Y ,1/2-Z \t",
        " ( 5) 1/2-Z ,1/2+X ,    Y \t( 6)     Y ,1/2-Z ,1/2+X \t",
        " ( 7)    -Z ,1/2+X ,1/2-Y \t( 8) 1/2-Y ,   -Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/2-Z ,   -X \t(10)    -X ,1/2+Y ,1/2-Z \t",
        " (11) 1/2-Z ,   -X ,1/2+Y \t(12) 1/2+X ,1/2-Y ,   -Z \t",
        " (13) 1/4+Y ,1/4+X ,1/4+Z \t(14) 1/4+Z ,1/4+Y ,1/4+X \t",
        " (15) 1/4+X ,1/4+Z ,1/4+Y \t(16) 3/4+Y ,1/4+X ,1/4-Z \t",
        " (17) 1/4-Z ,3/4+Y ,1/4+X \t(18) 1/4+X ,1/4-Z ,3/4+Y \t",
        " (19) 3/4-Z ,3/4+Y ,1/4-X \t(20) 1/4-X ,3/4-Z ,3/4+Y \t",
        " (21) 3/4+X ,1/4-Z ,3/4-Y \t(22) 3/4-Y ,3/4+X ,1/4-Z \t",
        " (23) 1/4-Z ,3/4-Y ,3/4+X \t(24) 3/4+Y ,1/4-X ,3/4-Z \t",
        " ",
    ],
    "p a -3": [
        " Space Group: P a -3",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,    Y ,1/2-Z \t",
        " ( 5) 1/2-Z ,1/2+X ,    Y \t( 6)     Y ,1/2-Z ,1/2+X \t",
        " ( 7)    -Z ,1/2+X ,1/2-Y \t( 8) 1/2-Y ,   -Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/2-Z ,   -X \t(10)    -X ,1/2+Y ,1/2-Z \t",
        " (11) 1/2-Z ,   -X ,1/2+Y \t(12) 1/2+X ,1/2-Y ,   -Z \t",
        " ",
    ],
    "c 2": [
        " Space Group: C 2",
        " The lattice is noncentrosymmetric C-centered monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in y",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "r 3 2 h": [
        " Space Group: R 3 2 h",
        " The lattice is noncentrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 18",
        " The Laue symmetry is 3m1",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,   Y-X,   -Z \t\t( 6)    X-Y,   -Y ,   -Z \t\t",
        " ",
    ],
    "p c": [
        " Space Group: P c",
        " The lattice is noncentrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 2",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in x z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "f 2 2 2": [
        " Space Group: F 2 2 2",
        " The lattice is noncentrosymmetric F-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ( 3)    -X ,    Y ,   -Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 6 c c": [
        " Space Group: P 6 c c",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ( 7)    Y-X,    Y ,1/2+Z \t( 8)    -X ,   Y-X,1/2+Z \t",
        " ( 9)    -Y ,   -X ,1/2+Z \t(10)    X-Y,   -Y ,1/2+Z \t",
        " (11)     X ,   X-Y,1/2+Z \t(12)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "r 3 2 r": [
        " Space Group: R 3 2 r",
        " The lattice is noncentrosymmetric primitive rhombohedral",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3mR",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)    -Y ,   -X ,   -Z \t\t",
        " ( 5)    -Z ,   -Y ,   -X \t\t( 6)    -X ,   -Z ,   -Y \t\t",
        " ",
    ],
    "p 62 2 2": [
        " Space Group: P 62 2 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/3+Z \t",
        " ( 3)    -Y ,   X-Y,2/3+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,1/3+Z \t( 6)     Y ,   Y-X,2/3+Z \t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)     X ,   X-Y,1/3-Z \t",
        " ( 9)     Y ,    X ,2/3-Z \t(10)    Y-X,    Y ,   -Z \t\t",
        " (11)    -X ,   Y-X,1/3-Z \t(12)    -Y ,   -X ,2/3-Z \t",
        " ",
    ],
    "i 41 c d": [
        " Space Group: I 41 c d",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,1/2+X ,1/4+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,1/2+Z \t( 4) 1/2+Y ,   -X ,3/4+Z \t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)    -Y ,1/2-X ,3/4+Z \t",
        " ( 7) 1/2+X ,1/2-Y ,    Z \t( 8) 1/2+Y ,    X ,1/4+Z \t",
        " ",
    ],
    "f m m 2": [
        " Space Group: F m m 2",
        " The lattice is noncentrosymmetric F-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p m m 2": [
        " Space Group: P m m 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p m": [
        " Space Group: P m",
        " The lattice is noncentrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 2",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in x z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,    Z \t\t",
        " ",
    ],
    "i 4 2 2": [
        " Space Group: I 4 2 2",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,   -Z \t\t( 6)    -Y ,   -X ,   -Z \t\t",
        " ( 7)     X ,   -Y ,   -Z \t\t( 8)     Y ,    X ,   -Z \t\t",
        " ",
    ],
    "p 31 2 1": [
        " Space Group: P 31 2 1",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3m1",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,1/3+Z \t",
        " ( 3)    Y-X,   -X ,2/3+Z \t( 4)     Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,   Y-X,1/3-Z \t( 6)    X-Y,   -Y ,2/3-Z \t",
        " ",
    ],
    "i -4": [
        " Space Group: I -4",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ",
    ],
    "p 2": [
        " Space Group: P 2",
        " The lattice is noncentrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 2",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in y",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "p 21 21 2": [
        " Space Group: P 21 21 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2+X ,1/2-Y ,   -Z \t",
        " ( 3) 1/2-X ,1/2+Y ,   -Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "c 1 2/c 1": [
        " Space Group: C 1 2/c 1",
        " The lattice is centrosymmetric C-centered monoclinic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,1/2-Z \t",
        " ",
    ],
    "i b a 2": [
        " Space Group: I b a 2",
        " The lattice is noncentrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3) 1/2+X ,1/2-Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p b a m": [
        " Space Group: P b a m",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3) 1/2+X ,1/2-Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p b a n": [
        " Space Group: P b a n",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,    Z \t",
        " ( 3) 1/2+X ,   -Y ,    Z \t( 4) 1/2-X ,1/2-Y ,    Z \t",
        " ",
    ],
    "p m -3 n": [
        " Space Group: P m -3 n",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " (13) 1/2+Y ,1/2+X ,1/2+Z \t(14) 1/2+Z ,1/2+Y ,1/2+X \t",
        " (15) 1/2+X ,1/2+Z ,1/2+Y \t(16) 1/2+Y ,1/2+X ,1/2-Z \t",
        " (17) 1/2-Z ,1/2+Y ,1/2+X \t(18) 1/2+X ,1/2-Z ,1/2+Y \t",
        " (19) 1/2-Z ,1/2+Y ,1/2-X \t(20) 1/2-X ,1/2-Z ,1/2+Y \t",
        " (21) 1/2+X ,1/2-Z ,1/2-Y \t(22) 1/2-Y ,1/2+X ,1/2-Z \t",
        " (23) 1/2-Z ,1/2-Y ,1/2+X \t(24) 1/2+Y ,1/2-X ,1/2-Z \t",
        " ",
    ],
    "i b c a": [
        " Space Group: I b c a",
        " The lattice is centrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3)     X ,1/2-Y ,1/2+Z \t( 4) 1/2-X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 42 21 2": [
        " Space Group: P 42 21 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,1/2+X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4) 1/2+Y ,1/2-X ,1/2+Z \t",
        " ( 5) 1/2-X ,1/2+Y ,1/2-Z \t( 6)    -Y ,   -X ,   -Z \t\t",
        " ( 7) 1/2+X ,1/2-Y ,1/2-Z \t( 8)     Y ,    X ,   -Z \t\t",
        " ",
    ],
    "f m m m": [
        " Space Group: F m m m",
        " The lattice is centrosymmetric F-centered orthorhombic",
        " Multiplicity of a general site is 32",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "i 41": [
        " Space Group: I 41",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,1/2+X ,1/4+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,1/2+Z \t( 4) 1/2+Y ,   -X ,3/4+Z \t",
        " ",
    ],
    "p 6": [
        " Space Group: P 6",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ",
    ],
    "p 3": [
        " Space Group: P 3",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 3",
        " The Laue symmetry is 3",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t",
    ],
    "p m m n": [
        " Space Group: P m m n",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,    Z \t",
        " ( 3)     X ,1/2-Y ,    Z \t( 4) 1/2-X ,1/2-Y ,    Z \t",
        " ",
    ],
    "p m m m": [
        " Space Group: P m m m",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "r 3 m h": [
        " Space Group: R 3 m h",
        " The lattice is noncentrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 18",
        " The Laue symmetry is 3m1",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,    Z \t\t",
        " ( 5)    -Y ,   -X ,    Z \t\t( 6)     X ,   X-Y,    Z \t\t",
        " ",
    ],
    "c m c 21": [
        " Space Group: C m c 21",
        " The lattice is noncentrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 1": [
        " Space Group: P 1",
        " The lattice is noncentrosymmetric primitive triclinic",
        " Multiplicity of a general site is 1",
        " The Laue symmetry is -1",
        " The location of the origin is arbitrary in xyz",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t",
    ],
    "i 4": [
        " Space Group: I 4",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ",
    ],
    "p 4": [
        " Space Group: P 4",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 4/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ",
    ],
    "p 42 b c": [
        " Space Group: P 42 b c",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ( 5) 1/2-X ,1/2+Y ,    Z \t( 6) 1/2-Y ,1/2-X ,1/2+Z \t",
        " ( 7) 1/2+X ,1/2-Y ,    Z \t( 8) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ",
    ],
    "p m m a": [
        " Space Group: P m m a",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,    Z \t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4) 1/2-X ,   -Y ,    Z \t",
        " ",
    ],
    "i -4 m 2": [
        " Space Group: I -4 m 2",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)     Y ,    X ,   -Z \t\t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)    -Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "p -4 21 c": [
        " Space Group: P -4 21 c",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,1/2-Z \t( 6) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ( 7) 1/2+X ,1/2-Y ,1/2-Z \t( 8) 1/2-Y ,1/2-X ,1/2+Z \t",
        " ",
    ],
    "p 4 2 2": [
        " Space Group: P 4 2 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,   -Z \t\t( 6)    -Y ,   -X ,   -Z \t\t",
        " ( 7)     X ,   -Y ,   -Z \t\t( 8)     Y ,    X ,   -Z \t\t",
        " ",
    ],
    "p 4/m c c": [
        " Space Group: P 4/m c c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p -6 2 m": [
        " Space Group: P -6 2 m",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    Y-X,   -X ,   -Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)    -Y ,   X-Y,   -Z \t\t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)    -X ,   Y-X,    Z \t\t",
        " ( 9)     Y ,    X ,   -Z \t\t(10)    X-Y,   -Y ,    Z \t\t",
        " (11)    -X ,   Y-X,   -Z \t\t(12)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p 1 2/m 1": [
        " Space Group: P 1 2/m 1",
        " The lattice is centrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "p -6 2 c": [
        " Space Group: P -6 2 c",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    Y-X,   -X ,1/2-Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)     X ,    Y ,1/2-Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)    -Y ,   X-Y,1/2-Z \t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)    -X ,   Y-X,1/2+Z \t",
        " ( 9)     Y ,    X ,   -Z \t\t(10)    X-Y,   -Y ,1/2+Z \t",
        " (11)    -X ,   Y-X,   -Z \t\t(12)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p 6 m m": [
        " Space Group: P 6 m m",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ( 7)    Y-X,    Y ,    Z \t\t( 8)    -X ,   Y-X,    Z \t\t",
        " ( 9)    -Y ,   -X ,    Z \t\t(10)    X-Y,   -Y ,    Z \t\t",
        " (11)     X ,   X-Y,    Z \t\t(12)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "c c": [
        " Space Group: C c",
        " The lattice is noncentrosymmetric C-centered monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in x z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 43 3 2": [
        " Space Group: P 43 3 2",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 3/4-Y ,1/4+X ,3/4+Z \t",
        " ( 5) 3/4+Z ,3/4-Y ,1/4+X \t( 6) 1/4+X ,3/4+Z ,3/4-Y \t",
        " ( 7) 1/2-X ,   -Y ,1/2+Z \t( 8)    -Z ,1/2+X ,1/2-Y \t",
        " ( 9) 1/2-Y ,   -Z ,1/2+X \t(10) 1/2+X ,1/2-Y ,   -Z \t",
        " (11) 1/2+Z ,1/2-X ,   -Y \t(12)    -Y ,1/2+Z ,1/2-X \t",
        " (13) 3/4+Y ,3/4-X ,1/4+Z \t(14) 1/4+Z ,3/4+Y ,3/4-X \t",
        " (15) 3/4-X ,1/4+Z ,3/4+Y \t(16) 1/4-X ,1/4-Z ,1/4-Y \t",
        " (17) 1/4-Y ,1/4-X ,1/4-Z \t(18) 1/4-Z ,1/4-Y ,1/4-X \t",
        " (19) 1/2+Y ,1/2-Z ,   -X \t(20) 1/4+Y ,3/4+X ,3/4-Z \t",
        " (21) 3/4-Z ,1/4+Y ,3/4+X \t(22) 3/4+X ,3/4-Z ,1/4+Y \t",
        " (23)    -X ,1/2+Y ,1/2-Z \t(24) 1/2-Z ,   -X ,1/2+Y \t",
        " ",
    ],
    "a 2 2 2": [
        " Space Group: A 2 2 2",
        " The lattice is noncentrosymmetric A-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ( 3)    -X ,    Y ,   -Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p -3": [
        " Space Group: P -3",
        " The lattice is centrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t",
    ],
    "p -1": [
        " Space Group: P -1",
        " The lattice is centrosymmetric primitive triclinic",
        " Multiplicity of a general site is 2",
        " The Laue symmetry is -1",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t",
    ],
    "f d d 2": [
        " Space Group: F d d 2",
        " The lattice is noncentrosymmetric F-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/4-X ,1/4+Y ,1/4+Z \t",
        " ( 3) 1/4+X ,1/4-Y ,1/4+Z \t( 4)    -X ,1/2-Y ,1/2+Z \t",
        " ",
    ],
    "p 62": [
        " Space Group: P 62",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/3+Z \t",
        " ( 3)    -Y ,   X-Y,2/3+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,1/3+Z \t( 6)     Y ,   Y-X,2/3+Z \t",
        " ",
    ],
    "c m m a": [
        " Space Group: C m m a",
        " The lattice is centrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,1/2-Y ,    Z \t( 4)    -X ,1/2-Y ,    Z \t",
        " ",
    ],
    "p -3 c 1": [
        " Space Group: P -3 c 1",
        " The lattice is centrosymmetric primitive trigonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 3m1",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,1/2+Z \t",
        " ( 5)    -Y ,   -X ,1/2+Z \t( 6)     X ,   X-Y,1/2+Z \t",
        " ",
    ],
    "c m c m": [
        " Space Group: C m c m",
        " The lattice is centrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "c m m m": [
        " Space Group: C m m m",
        " The lattice is centrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "c m c a": [
        " Space Group: C m c a",
        " The lattice is centrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,1/2-Y ,1/2+Z \t( 4)    -X ,1/2-Y ,1/2+Z \t",
        " ",
    ],
    "i a -3": [
        " Space Group: I a -3",
        " The lattice is centrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,    Y ,1/2-Z \t",
        " ( 5) 1/2-Z ,1/2+X ,    Y \t( 6)     Y ,1/2-Z ,1/2+X \t",
        " ( 7)    -Z ,1/2+X ,1/2-Y \t( 8) 1/2-Y ,   -Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/2-Z ,   -X \t(10)    -X ,1/2+Y ,1/2-Z \t",
        " (11) 1/2-Z ,   -X ,1/2+Y \t(12) 1/2+X ,1/2-Y ,   -Z \t",
        " ",
    ],
    "i m a 2": [
        " Space Group: I m a 2",
        " The lattice is noncentrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,    Z \t",
        " ( 3) 1/2+X ,   -Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 63/m m c": [
        " Space Group: P 63/m m c",
        " The lattice is centrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is 6/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/2+Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,1/2+Z \t",
        " ( 7)    Y-X,    Y ,    Z \t\t( 8)    -X ,   Y-X,1/2+Z \t",
        " ( 9)    -Y ,   -X ,    Z \t\t(10)    X-Y,   -Y ,1/2+Z \t",
        " (11)     X ,   X-Y,    Z \t\t(12)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p 4 3 2": [
        " Space Group: P 4 3 2",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)    -Y ,    X ,    Z \t\t",
        " ( 5)     Z ,   -Y ,    X \t\t( 6)     X ,    Z ,   -Y \t\t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13)     Y ,   -X ,    Z \t\t(14)     Z ,    Y ,   -X \t\t",
        " (15)    -X ,    Z ,    Y \t\t(16)    -X ,   -Z ,   -Y \t\t",
        " (17)    -Y ,   -X ,   -Z \t\t(18)    -Z ,   -Y ,   -X \t\t",
        " (19)     Y ,   -Z ,   -X \t\t(20)     Y ,    X ,   -Z \t\t",
        " (21)    -Z ,    Y ,    X \t\t(22)     X ,   -Z ,    Y \t\t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "i 41/a c d": [
        " Space Group: I 41/a c d",
        " The lattice is centrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 32",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/4-Y ,3/4+X ,1/4+Z \t",
        " ( 3) 1/2-X ,   -Y ,1/2+Z \t( 4) 1/4+Y ,1/4-X ,3/4+Z \t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6) 1/4-Y ,3/4-X ,3/4+Z \t",
        " ( 7) 1/2+X ,   -Y ,    Z \t( 8) 1/4+Y ,1/4+X ,1/4+Z \t",
        " ",
    ],
    "p -4 2 c": [
        " Space Group: P -4 2 c",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,    Y ,1/2-Z \t( 6)     Y ,    X ,1/2+Z \t",
        " ( 7)     X ,   -Y ,1/2-Z \t( 8)    -Y ,   -X ,1/2+Z \t",
        " ",
    ],
    "p 4 n c": [
        " Space Group: P 4 n c",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,1/2+Z \t( 6) 1/2-Y ,1/2-X ,1/2+Z \t",
        " ( 7) 1/2+X ,1/2-Y ,1/2+Z \t( 8) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ",
    ],
    "p 4/m": [
        " Space Group: P 4/m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ",
    ],
    "p 4/n": [
        " Space Group: P 4/n",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,    Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,    Z \t",
        " ",
    ],
    "p 21/c": [
        " Space Group: P 21/c",
        " The lattice is centrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,1/2-Z \t",
        " ",
    ],
    "r -3 c": [
        " Space Group: R -3 c",
        " The lattice is centrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 36",
        " The Laue symmetry is 3m1",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,1/2+Z \t",
        " ( 5)    -Y ,   -X ,1/2+Z \t( 6)     X ,   X-Y,1/2+Z \t",
        " ",
    ],
    "p 4/n m m": [
        " Space Group: P 4/n m m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,    Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,    Z \t",
        " ( 5) 1/2-X ,    Y ,    Z \t( 6) 1/2-Y ,1/2-X ,    Z \t",
        " ( 7)     X ,1/2-Y ,    Z \t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "i 4/m": [
        " Space Group: I 4/m",
        " The lattice is centrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ",
    ],
    "p 3 m 1": [
        " Space Group: P 3 m 1",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3m1",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,    Z \t\t",
        " ( 5)    -Y ,   -X ,    Z \t\t( 6)     X ,   X-Y,    Z \t\t",
        " ",
    ],
    "p 63/m": [
        " Space Group: P 63/m",
        " The lattice is centrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/2+Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,1/2+Z \t",
        " ",
    ],
    "p 6 2 2": [
        " Space Group: P 6 2 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)     X ,   X-Y,   -Z \t\t",
        " ( 9)     Y ,    X ,   -Z \t\t(10)    Y-X,    Y ,   -Z \t\t",
        " (11)    -X ,   Y-X,   -Z \t\t(12)    -Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "p 2/m": [
        " Space Group: P 2/m",
        " The lattice is centrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "f d 3": [
        " Space Group: F d 3",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/4+X ,1/4+Y ,   -Z \t",
        " ( 5)    -Z ,1/4+X ,1/4+Y \t( 6) 1/4+Y ,   -Z ,1/4+X \t",
        " ( 7) 1/4-Z ,1/2+X ,3/4-Y \t( 8) 3/4-Y ,1/4-Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/4-Z ,3/4-X \t(10) 3/4-X ,1/2+Y ,1/4-Z \t",
        " (11) 1/4-Z ,3/4-X ,1/2+Y \t(12) 1/2+X ,3/4-Y ,1/4-Z \t",
        " ",
    ],
    "i 41/a m d": [
        " Space Group: I 41/a m d",
        " The lattice is centrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 32",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/4-Y ,3/4+X ,1/4+Z \t",
        " ( 3) 1/2-X ,   -Y ,1/2+Z \t( 4) 1/4+Y ,1/4-X ,3/4+Z \t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6) 1/4-Y ,3/4-X ,1/4+Z \t",
        " ( 7) 1/2+X ,   -Y ,1/2+Z \t( 8) 1/4+Y ,1/4+X ,3/4+Z \t",
        " ",
    ],
    "p 4/n c c": [
        " Space Group: P 4/n c c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,    Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,    Z \t",
        " ( 5) 1/2-X ,    Y ,1/2+Z \t( 6) 1/2-Y ,1/2-X ,1/2+Z \t",
        " ( 7)     X ,1/2-Y ,1/2+Z \t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "i m m a": [
        " Space Group: I m m a",
        " The lattice is centrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,1/2-Y ,    Z \t( 4)    -X ,1/2-Y ,    Z \t",
        " ",
    ],
    "p 4 b m": [
        " Space Group: P 4 b m",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,    Z \t( 6) 1/2-Y ,1/2-X ,    Z \t",
        " ( 7) 1/2+X ,1/2-Y ,    Z \t( 8) 1/2+Y ,1/2+X ,    Z \t",
        " ",
    ],
    "p 2/c": [
        " Space Group: P 2/c",
        " The lattice is centrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,1/2-Z \t",
        " ",
    ],
    "p -6 m 2": [
        " Space Group: P -6 m 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    Y-X,   -X ,   -Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)    -Y ,   X-Y,   -Z \t\t",
        " ( 7)    Y-X,    Y ,    Z \t\t( 8)     X ,   X-Y,   -Z \t\t",
        " ( 9)    -Y ,   -X ,    Z \t\t(10)    Y-X,    Y ,   -Z \t\t",
        " (11)     X ,   X-Y,    Z \t\t(12)    -Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "p n n 2": [
        " Space Group: P n n 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,1/2+Z \t",
        " ( 3) 1/2+X ,1/2-Y ,1/2+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 31 1 2": [
        " Space Group: P 31 1 2",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 31m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,1/3+Z \t",
        " ( 3)    Y-X,   -X ,2/3+Z \t( 4)     X ,   X-Y,   -Z \t\t",
        " ( 5)    Y-X,    Y ,1/3-Z \t( 6)    -Y ,   -X ,2/3-Z \t",
        " ",
    ],
    "f -4 3 c": [
        " Space Group: F -4 3 c",
        " The lattice is noncentrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+Y ,1/2-X ,1/2-Z \t",
        " ( 5) 1/2-Z ,1/2+Y ,1/2-X \t( 6) 1/2-X ,1/2-Z ,1/2+Y \t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13) 1/2-Y ,1/2+X ,1/2-Z \t(14) 1/2-Z ,1/2-Y ,1/2+X \t",
        " (15) 1/2+X ,1/2-Z ,1/2-Y \t(16) 1/2+X ,1/2+Z ,1/2+Y \t",
        " (17) 1/2+Y ,1/2+X ,1/2+Z \t(18) 1/2+Z ,1/2+Y ,1/2+X \t",
        " (19)     Y ,   -Z ,   -X \t\t(20) 1/2-Y ,1/2-X ,1/2+Z \t",
        " (21) 1/2+Z ,1/2-Y ,1/2-X \t(22) 1/2-X ,1/2+Z ,1/2-Y \t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "i m -3": [
        " Space Group: I m -3",
        " The lattice is centrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " ",
    ],
    "f -4 3 m": [
        " Space Group: F -4 3 m",
        " The lattice is noncentrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     Y ,   -X ,   -Z \t\t",
        " ( 5)    -Z ,    Y ,   -X \t\t( 6)    -X ,   -Z ,    Y \t\t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13)    -Y ,    X ,   -Z \t\t(14)    -Z ,   -Y ,    X \t\t",
        " (15)     X ,   -Z ,   -Y \t\t(16)     X ,    Z ,    Y \t\t",
        " (17)     Y ,    X ,    Z \t\t(18)     Z ,    Y ,    X \t\t",
        " (19)     Y ,   -Z ,   -X \t\t(20)    -Y ,   -X ,    Z \t\t",
        " (21)     Z ,   -Y ,   -X \t\t(22)    -X ,    Z ,   -Y \t\t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "i 21 3": [
        " Space Group: I 21 3",
        " The lattice is noncentrosymmetric I-centered cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,1/2-Y ,   -Z \t",
        " ( 5)    -Z ,1/2+X ,1/2-Y \t( 6) 1/2-Y ,   -Z ,1/2+X \t",
        " ( 7) 1/2-Z ,   -X ,1/2+Y \t( 8) 1/2+Y ,1/2-Z ,   -X \t",
        " ( 9)    -Y ,1/2+Z ,1/2-X \t(10) 1/2-X ,   -Y ,1/2+Z \t",
        " (11) 1/2+Z ,1/2-X ,   -Y \t(12)    -X ,1/2+Y ,1/2-Z \t",
        " ",
    ],
    "p 42/m m c": [
        " Space Group: P 42/m m c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p 65 2 2": [
        " Space Group: P 65 2 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,5/6+Z \t",
        " ( 3)    -Y ,   X-Y,2/3+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,1/3+Z \t( 6)     Y ,   Y-X,1/6+Z \t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)     X ,   X-Y,5/6-Z \t",
        " ( 9)     Y ,    X ,2/3-Z \t(10)    Y-X,    Y ,1/2-Z \t",
        " (11)    -X ,   Y-X,1/3-Z \t(12)    -Y ,   -X ,1/6-Z \t",
        " ",
    ],
    "p 4/m n c": [
        " Space Group: P 4/m n c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,1/2+Z \t( 6) 1/2-Y ,1/2-X ,1/2+Z \t",
        " ( 7) 1/2+X ,1/2-Y ,1/2+Z \t( 8) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ",
    ],
    "c 2/m": [
        " Space Group: C 2/m",
        " The lattice is centrosymmetric C-centered monoclinic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "f d d d": [
        " Space Group: F d d d",
        " The lattice is centrosymmetric F-centered orthorhombic",
        " Multiplicity of a general site is 32",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/4+Y ,1/4+Z \t",
        " ( 3) 1/4+X ,   -Y ,1/4+Z \t( 4) 3/4-X ,1/4-Y ,1/2+Z \t",
        " ",
    ],
    "c m m 2": [
        " Space Group: C m m 2",
        " The lattice is noncentrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 43 21 2": [
        " Space Group: P 43 21 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,1/2+X ,3/4+Z \t",
        " ( 3)    -X ,   -Y ,1/2+Z \t( 4) 1/2+Y ,1/2-X ,1/4+Z \t",
        " ( 5) 1/2-X ,1/2+Y ,3/4-Z \t( 6)    -Y ,   -X ,1/2-Z \t",
        " ( 7) 1/2+X ,1/2-Y ,1/4-Z \t( 8)     Y ,    X ,   -Z \t\t",
        " ",
    ],
    "p -3 1 m": [
        " Space Group: P -3 1 m",
        " The lattice is centrosymmetric primitive trigonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 31m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     Y ,    X ,    Z \t\t",
        " ( 5)    -X ,   Y-X,    Z \t\t( 6)    X-Y,   -Y ,    Z \t\t",
        " ",
    ],
    "i 2 2 2": [
        " Space Group: I 2 2 2",
        " The lattice is noncentrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ( 3)    -X ,    Y ,   -Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 42/n b c": [
        " Space Group: P 42/n b c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,1/2+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,1/2+Z \t",
        " ( 5)    -X ,1/2+Y ,    Z \t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7) 1/2+X ,   -Y ,    Z \t( 8) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ",
    ],
    "i 4 3 2": [
        " Space Group: I 4 3 2",
        " The lattice is noncentrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)    -Y ,    X ,    Z \t\t",
        " ( 5)     Z ,   -Y ,    X \t\t( 6)     X ,    Z ,   -Y \t\t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13)     Y ,   -X ,    Z \t\t(14)     Z ,    Y ,   -X \t\t",
        " (15)    -X ,    Z ,    Y \t\t(16)    -X ,   -Z ,   -Y \t\t",
        " (17)    -Y ,   -X ,   -Z \t\t(18)    -Z ,   -Y ,   -X \t\t",
        " (19)     Y ,   -Z ,   -X \t\t(20)     Y ,    X ,   -Z \t\t",
        " (21)    -Z ,    Y ,    X \t\t(22)     X ,   -Z ,    Y \t\t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "p 41 3 2": [
        " Space Group: P 41 3 2",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/4-Y ,3/4+X ,1/4+Z \t",
        " ( 5) 1/4+Z ,1/4-Y ,3/4+X \t( 6) 3/4+X ,1/4+Z ,1/4-Y \t",
        " ( 7) 1/2-X ,   -Y ,1/2+Z \t( 8)    -Z ,1/2+X ,1/2-Y \t",
        " ( 9) 1/2-Y ,   -Z ,1/2+X \t(10) 1/2+X ,1/2-Y ,   -Z \t",
        " (11) 1/2+Z ,1/2-X ,   -Y \t(12)    -Y ,1/2+Z ,1/2-X \t",
        " (13) 1/4+Y ,1/4-X ,3/4+Z \t(14) 3/4+Z ,1/4+Y ,1/4-X \t",
        " (15) 1/4-X ,3/4+Z ,1/4+Y \t(16) 3/4-X ,3/4-Z ,3/4-Y \t",
        " (17) 3/4-Y ,3/4-X ,3/4-Z \t(18) 3/4-Z ,3/4-Y ,3/4-X \t",
        " (19) 1/2+Y ,1/2-Z ,   -X \t(20) 3/4+Y ,1/4+X ,1/4-Z \t",
        " (21) 1/4-Z ,3/4+Y ,1/4+X \t(22) 1/4+X ,1/4-Z ,3/4+Y \t",
        " (23)    -X ,1/2+Y ,1/2-Z \t(24) 1/2-Z ,   -X ,1/2+Y \t",
        " ",
    ],
    "p 42/n m c": [
        " Space Group: P 42/n m c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,1/2+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,1/2+Z \t",
        " ( 5) 1/2-X ,    Y ,    Z \t( 6) 1/2-Y ,1/2-X ,1/2+Z \t",
        " ( 7)     X ,1/2-Y ,    Z \t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p 64 2 2": [
        " Space Group: P 64 2 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,2/3+Z \t",
        " ( 3)    -Y ,   X-Y,1/3+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,2/3+Z \t( 6)     Y ,   Y-X,1/3+Z \t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)     X ,   X-Y,2/3-Z \t",
        " ( 9)     Y ,    X ,1/3-Z \t(10)    Y-X,    Y ,   -Z \t\t",
        " (11)    -X ,   Y-X,2/3-Z \t(12)    -Y ,   -X ,1/3-Z \t",
        " ",
    ],
    "p c a 21": [
        " Space Group: P c a 21",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,1/2+Z \t",
        " ( 3) 1/2+X ,   -Y ,    Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "f d -3 c": [
        " Space Group: F d -3 c",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 192",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/4+X ,1/4+Y ,   -Z \t",
        " ( 5)    -Z ,1/4+X ,1/4+Y \t( 6) 1/4+Y ,   -Z ,1/4+X \t",
        " ( 7) 1/4-Z ,1/2+X ,3/4-Y \t( 8) 3/4-Y ,1/4-Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/4-Z ,3/4-X \t(10) 3/4-X ,1/2+Y ,1/4-Z \t",
        " (11) 1/4-Z ,3/4-X ,1/2+Y \t(12) 1/2+X ,3/4-Y ,1/4-Z \t",
        " (13)     Y ,    X ,1/2+Z \t(14) 1/2+Z ,    Y ,    X \t",
        " (15)     X ,1/2+Z ,    Y \t(16) 1/4+Y ,1/4+X ,1/2-Z \t",
        " (17) 1/2-Z ,1/4+Y ,1/4+X \t(18) 1/4+X ,1/2-Z ,1/4+Y \t",
        " (19) 3/4-Z ,1/2+Y ,3/4-X \t(20) 3/4-X ,3/4-Z ,1/2+Y \t",
        " (21) 1/2+X ,3/4-Z ,3/4-Y \t(22) 3/4-Y ,1/2+X ,3/4-Z \t",
        " (23) 3/4-Z ,3/4-Y ,1/2+X \t(24) 1/2+Y ,3/4-X ,3/4-Z \t",
        " ",
    ],
    "p n a 21": [
        " Space Group: P n a 21",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,1/2+Z \t",
        " ( 3) 1/2+X ,1/2-Y ,    Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p -4 n 2": [
        " Space Group: P -4 n 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,1/2+Z \t( 6) 1/2+Y ,1/2+X ,1/2-Z \t",
        " ( 7) 1/2+X ,1/2-Y ,1/2+Z \t( 8) 1/2-Y ,1/2-X ,1/2-Z \t",
        " ",
    ],
    "p 42/n n m": [
        " Space Group: P 42/n n m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,1/2+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,1/2+Z \t",
        " ( 5)    -X ,1/2+Y ,1/2+Z \t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7) 1/2+X ,   -Y ,1/2+Z \t( 8) 1/2+Y ,1/2+X ,    Z \t",
        " ",
    ],
    "f d -3 m": [
        " Space Group: F d -3 m",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 192",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/4+X ,1/4+Y ,   -Z \t",
        " ( 5)    -Z ,1/4+X ,1/4+Y \t( 6) 1/4+Y ,   -Z ,1/4+X \t",
        " ( 7) 1/4-Z ,1/2+X ,3/4-Y \t( 8) 3/4-Y ,1/4-Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/4-Z ,3/4-X \t(10) 3/4-X ,1/2+Y ,1/4-Z \t",
        " (11) 1/4-Z ,3/4-X ,1/2+Y \t(12) 1/2+X ,3/4-Y ,1/4-Z \t",
        " (13)     Y ,    X ,    Z \t\t(14)     Z ,    Y ,    X \t\t",
        " (15)     X ,    Z ,    Y \t\t(16) 1/4+Y ,1/4+X ,   -Z \t",
        " (17)    -Z ,1/4+Y ,1/4+X \t(18) 1/4+X ,   -Z ,1/4+Y \t",
        " (19) 1/4-Z ,1/2+Y ,3/4-X \t(20) 3/4-X ,1/4-Z ,1/2+Y \t",
        " (21) 1/2+X ,1/4-Z ,3/4-Y \t(22) 3/4-Y ,1/2+X ,1/4-Z \t",
        " (23) 1/4-Z ,3/4-Y ,1/2+X \t(24) 1/2+Y ,3/4-X ,1/4-Z \t",
        " ",
    ],
    "r -3 c r": [
        " Space Group: R -3 c r",
        " The lattice is centrosymmetric primitive rhombohedral",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 3mR",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ( 5) 1/2+Z ,1/2+Y ,1/2+X \t( 6) 1/2+X ,1/2+Z ,1/2+Y \t",
        " ",
    ],
    "p 63 m c": [
        " Space Group: P 63 m c",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/2+Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,1/2+Z \t",
        " ( 7)    Y-X,    Y ,    Z \t\t( 8)    -X ,   Y-X,1/2+Z \t",
        " ( 9)    -Y ,   -X ,    Z \t\t(10)    X-Y,   -Y ,1/2+Z \t",
        " (11)     X ,   X-Y,    Z \t\t(12)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p 4/m b m": [
        " Space Group: P 4/m b m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,    Z \t( 6) 1/2-Y ,1/2-X ,    Z \t",
        " ( 7) 1/2+X ,1/2-Y ,    Z \t( 8) 1/2+Y ,1/2+X ,    Z \t",
        " ",
    ],
    "p 2 2 2": [
        " Space Group: P 2 2 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ( 3)    -X ,    Y ,   -Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 63 2 2": [
        " Space Group: P 63 2 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/2+Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,1/2+Z \t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)     X ,   X-Y,1/2-Z \t",
        " ( 9)     Y ,    X ,   -Z \t\t(10)    Y-X,    Y ,1/2-Z \t",
        " (11)    -X ,   Y-X,   -Z \t\t(12)    -Y ,   -X ,1/2-Z \t",
        " ",
    ],
    "p 6/m m m": [
        " Space Group: P 6/m m m",
        " The lattice is centrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is 6/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ( 7)    Y-X,    Y ,    Z \t\t( 8)    -X ,   Y-X,    Z \t\t",
        " ( 9)    -Y ,   -X ,    Z \t\t(10)    X-Y,   -Y ,    Z \t\t",
        " (11)     X ,   X-Y,    Z \t\t(12)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p c c n": [
        " Space Group: P c c n",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,1/2+Z \t",
        " ( 3)     X ,1/2-Y ,1/2+Z \t( 4) 1/2-X ,1/2-Y ,    Z \t",
        " ",
    ],
    "p c c m": [
        " Space Group: P c c m",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,1/2+Z \t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p m n a": [
        " Space Group: P m n a",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3) 1/2+X ,   -Y ,1/2+Z \t( 4) 1/2-X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "f 41 3 2": [
        " Space Group: F 41 3 2",
        " The lattice is noncentrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 3/4-Y ,3/4+X ,1/4+Z \t",
        " ( 5) 1/4+Z ,3/4-Y ,3/4+X \t( 6) 3/4+X ,1/4+Z ,3/4-Y \t",
        " ( 7)    -X ,1/2-Y ,1/2+Z \t( 8) 1/2-Z ,1/2+X ,   -Y \t",
        " ( 9)    -Y ,1/2-Z ,1/2+X \t(10) 1/2+X ,   -Y ,1/2-Z \t",
        " (11) 1/2+Z ,   -X ,1/2-Y \t(12) 1/2-Y ,1/2+Z ,   -X \t",
        " (13) 1/4+Y ,3/4-X ,3/4+Z \t(14) 3/4+Z ,1/4+Y ,3/4-X \t",
        " (15) 3/4-X ,3/4+Z ,1/4+Y \t(16) 1/4-X ,1/4-Z ,1/4-Y \t",
        " (17) 1/4-Y ,1/4-X ,1/4-Z \t(18) 1/4-Z ,1/4-Y ,1/4-X \t",
        " (19) 1/2+Y ,   -Z ,1/2-X \t(20) 3/4+Y ,1/4+X ,3/4-Z \t",
        " (21) 3/4-Z ,3/4+Y ,1/4+X \t(22) 1/4+X ,3/4-Z ,3/4+Y \t",
        " (23) 1/2-X ,1/2+Y ,   -Z \t(24)    -Z ,1/2-X ,1/2+Y \t",
        " ",
    ],
    "r -3 r": [
        " Space Group: R -3 r",
        " The lattice is centrosymmetric primitive rhombohedral",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3R",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t",
    ],
    "p 1 1 2/m": [
        " Space Group: P 1 1 2/m",
        " The lattice is centrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is c",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 64": [
        " Space Group: P 64",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,2/3+Z \t",
        " ( 3)    -Y ,   X-Y,1/3+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,2/3+Z \t( 6)     Y ,   Y-X,1/3+Z \t",
        " ",
    ],
    "p c c a": [
        " Space Group: P c c a",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,1/2+Z \t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4) 1/2-X ,   -Y ,    Z \t",
        " ",
    ],
    "f m -3": [
        " Space Group: F m -3",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " ",
    ],
    "p -6": [
        " Space Group: P -6",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    Y-X,   -X ,   -Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)    -Y ,   X-Y,   -Z \t\t",
        " ",
    ],
    "i m m m": [
        " Space Group: I m m m",
        " The lattice is centrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p -4 2 m": [
        " Space Group: P -4 2 m",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,    Y ,   -Z \t\t( 6)     Y ,    X ,    Z \t\t",
        " ( 7)     X ,   -Y ,   -Z \t\t( 8)    -Y ,   -X ,    Z \t\t",
        " ",
    ],
    "p 21 3": [
        " Space Group: P 21 3",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is m3",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,1/2-Y ,   -Z \t",
        " ( 5)    -Z ,1/2+X ,1/2-Y \t( 6) 1/2-Y ,   -Z ,1/2+X \t",
        " ( 7) 1/2-Z ,   -X ,1/2+Y \t( 8) 1/2+Y ,1/2-Z ,   -X \t",
        " ( 9)    -Y ,1/2+Z ,1/2-X \t(10) 1/2-X ,   -Y ,1/2+Z \t",
        " (11) 1/2+Z ,1/2-X ,   -Y \t(12)    -X ,1/2+Y ,1/2-Z \t",
        " ",
    ],
    "p 4 m m": [
        " Space Group: P 4 m m",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p -4 m 2": [
        " Space Group: P -4 m 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)     Y ,    X ,   -Z \t\t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)    -Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "c 2/c": [
        " Space Group: C 2/c",
        " The lattice is centrosymmetric C-centered monoclinic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,1/2-Z \t",
        " ",
    ],
    "p 42 3 2": [
        " Space Group: P 42 3 2",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2-Y ,1/2+X ,1/2+Z \t",
        " ( 5) 1/2+Z ,1/2-Y ,1/2+X \t( 6) 1/2+X ,1/2+Z ,1/2-Y \t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13) 1/2+Y ,1/2-X ,1/2+Z \t(14) 1/2+Z ,1/2+Y ,1/2-X \t",
        " (15) 1/2-X ,1/2+Z ,1/2+Y \t(16) 1/2-X ,1/2-Z ,1/2-Y \t",
        " (17) 1/2-Y ,1/2-X ,1/2-Z \t(18) 1/2-Z ,1/2-Y ,1/2-X \t",
        " (19)     Y ,   -Z ,   -X \t\t(20) 1/2+Y ,1/2+X ,1/2-Z \t",
        " (21) 1/2-Z ,1/2+Y ,1/2+X \t(22) 1/2+X ,1/2-Z ,1/2+Y \t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "p 6/m c c": [
        " Space Group: P 6/m c c",
        " The lattice is centrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is 6/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ( 7)    Y-X,    Y ,1/2+Z \t( 8)    -X ,   Y-X,1/2+Z \t",
        " ( 9)    -Y ,   -X ,1/2+Z \t(10)    X-Y,   -Y ,1/2+Z \t",
        " (11)     X ,   X-Y,1/2+Z \t(12)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "f m 3": [
        " Space Group: F m 3",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " ",
    ],
    "p n n a": [
        " Space Group: P n n a",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,1/2+Z \t",
        " ( 3) 1/2+X ,1/2-Y ,1/2+Z \t( 4) 1/2-X ,   -Y ,    Z \t",
        " ",
    ],
    "i -4 3 d": [
        " Space Group: I -4 3 d",
        " The lattice is noncentrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 3/4+Y ,1/4-X ,3/4-Z \t",
        " ( 5) 3/4-Z ,3/4+Y ,1/4-X \t( 6) 1/4-X ,3/4-Z ,3/4+Y \t",
        " ( 7)    -X ,1/2-Y ,    Z \t( 8) 1/2-Z ,    X ,   -Y \t",
        " ( 9)    -Y ,1/2-Z ,    X \t(10)     X ,   -Y ,1/2-Z \t",
        " (11)     Z ,   -X ,1/2-Y \t(12) 1/2-Y ,    Z ,   -X \t",
        " (13) 1/4-Y ,1/4+X ,3/4-Z \t(14) 3/4-Z ,1/4-Y ,1/4+X \t",
        " (15) 1/4+X ,3/4-Z ,1/4-Y \t(16) 3/4+X ,3/4+Z ,3/4+Y \t",
        " (17) 3/4+Y ,3/4+X ,3/4+Z \t(18) 3/4+Z ,3/4+Y ,3/4+X \t",
        " (19) 1/2+Y ,1/2-Z ,   -X \t(20) 3/4-Y ,1/4-X ,1/4+Z \t",
        " (21) 1/4+Z ,3/4-Y ,1/4-X \t(22) 1/4-X ,1/4+Z ,3/4-Y \t",
        " (23)    -X ,1/2+Y ,1/2-Z \t(24) 1/2-Z ,   -X ,1/2+Y \t",
        " ",
    ],
    "p n n n": [
        " Space Group: P n n n",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,1/2+Z \t",
        " ( 3) 1/2+X ,   -Y ,1/2+Z \t( 4) 1/2-X ,1/2-Y ,    Z \t",
        " ",
    ],
    "p n n m": [
        " Space Group: P n n m",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,1/2+Z \t",
        " ( 3) 1/2+X ,1/2-Y ,1/2+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p -4": [
        " Space Group: P -4",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 4/m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ",
    ],
    "i -4 3 m": [
        " Space Group: I -4 3 m",
        " The lattice is noncentrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     Y ,   -X ,   -Z \t\t",
        " ( 5)    -Z ,    Y ,   -X \t\t( 6)    -X ,   -Z ,    Y \t\t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13)    -Y ,    X ,   -Z \t\t(14)    -Z ,   -Y ,    X \t\t",
        " (15)     X ,   -Z ,   -Y \t\t(16)     X ,    Z ,    Y \t\t",
        " (17)     Y ,    X ,    Z \t\t(18)     Z ,    Y ,    X \t\t",
        " (19)     Y ,   -Z ,   -X \t\t(20)    -Y ,   -X ,    Z \t\t",
        " (21)     Z ,   -Y ,   -X \t\t(22)    -X ,    Z ,   -Y \t\t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "p 65": [
        " Space Group: P 65",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,5/6+Z \t",
        " ( 3)    -Y ,   X-Y,2/3+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,1/3+Z \t( 6)     Y ,   Y-X,1/6+Z \t",
        " ",
    ],
    "r 3 r": [
        " Space Group: R 3 r",
        " The lattice is noncentrosymmetric primitive rhombohedral",
        " Multiplicity of a general site is 3",
        " The Laue symmetry is 3R",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t",
    ],
    "p 2/m 1 1": [
        " Space Group: P 2/m 1 1",
        " The lattice is centrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is a",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ",
    ],
    "i 41/a": [
        " Space Group: I 41/a",
        " The lattice is centrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 3/4-Y ,1/4+X ,1/4+Z \t",
        " ( 3) 1/2-X ,   -Y ,1/2+Z \t( 4) 3/4+Y ,3/4-X ,3/4+Z \t",
        " ",
    ],
    "p 63 c m": [
        " Space Group: P 63 c m",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/2+Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,1/2+Z \t",
        " ( 7)    Y-X,    Y ,1/2+Z \t( 8)    -X ,   Y-X,    Z \t\t",
        " ( 9)    -Y ,   -X ,1/2+Z \t(10)    X-Y,   -Y ,    Z \t\t",
        " (11)     X ,   X-Y,1/2+Z \t(12)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "c 1 2 1": [
        " Space Group: C 1 2 1",
        " The lattice is noncentrosymmetric C-centered monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in y",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "p b c n": [
        " Space Group: P b c n",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4) 1/2-X ,1/2-Y ,1/2+Z \t",
        " ",
    ],
    "p b c m": [
        " Space Group: P b c m",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,    Z \t",
        " ( 3)     X ,1/2-Y ,1/2+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "a m m 2": [
        " Space Group: A m m 2",
        " The lattice is noncentrosymmetric A-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "i m -3 m": [
        " Space Group: I m -3 m",
        " The lattice is centrosymmetric I-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " (13)     Y ,    X ,    Z \t\t(14)     Z ,    Y ,    X \t\t",
        " (15)     X ,    Z ,    Y \t\t(16)     Y ,    X ,   -Z \t\t",
        " (17)    -Z ,    Y ,    X \t\t(18)     X ,   -Z ,    Y \t\t",
        " (19)    -Z ,    Y ,   -X \t\t(20)    -X ,   -Z ,    Y \t\t",
        " (21)     X ,   -Z ,   -Y \t\t(22)    -Y ,    X ,   -Z \t\t",
        " (23)    -Z ,   -Y ,    X \t\t(24)     Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "i 4 m m": [
        " Space Group: I 4 m m",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p 61 2 2": [
        " Space Group: P 61 2 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/6+Z \t",
        " ( 3)    -Y ,   X-Y,1/3+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,2/3+Z \t( 6)     Y ,   Y-X,5/6+Z \t",
        " ( 7)    X-Y,   -Y ,   -Z \t\t( 8)     X ,   X-Y,1/6-Z \t",
        " ( 9)     Y ,    X ,1/3-Z \t(10)    Y-X,    Y ,1/2-Z \t",
        " (11)    -X ,   Y-X,2/3-Z \t(12)    -Y ,   -X ,5/6-Z \t",
        " ",
    ],
    "i m m 2": [
        " Space Group: I m m 2",
        " The lattice is noncentrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 42/n c m": [
        " Space Group: P 42/n c m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,1/2+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,1/2+Z \t",
        " ( 5) 1/2-X ,    Y ,1/2+Z \t( 6) 1/2-Y ,1/2-X ,    Z \t",
        " ( 7)     X ,1/2-Y ,1/2+Z \t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p b c a": [
        " Space Group: P b c a",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,    Z \t",
        " ( 3)     X ,1/2-Y ,1/2+Z \t( 4) 1/2-X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 4 21 2": [
        " Space Group: P 4 21 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,1/2+X ,    Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4) 1/2+Y ,1/2-X ,    Z \t",
        " ( 5) 1/2-X ,1/2+Y ,   -Z \t( 6)    -Y ,   -X ,   -Z \t\t",
        " ( 7) 1/2+X ,1/2-Y ,   -Z \t( 8)     Y ,    X ,   -Z \t\t",
        " ",
    ],
    "p 4/n n c": [
        " Space Group: P 4/n n c",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,    X ,    Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4)     Y ,1/2-X ,    Z \t",
        " ( 5)    -X ,1/2+Y ,1/2+Z \t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7) 1/2+X ,   -Y ,1/2+Z \t( 8) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ",
    ],
    "f m -3 m": [
        " Space Group: F m -3 m",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 192",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " (13)     Y ,    X ,    Z \t\t(14)     Z ,    Y ,    X \t\t",
        " (15)     X ,    Z ,    Y \t\t(16)     Y ,    X ,   -Z \t\t",
        " (17)    -Z ,    Y ,    X \t\t(18)     X ,   -Z ,    Y \t\t",
        " (19)    -Z ,    Y ,   -X \t\t(20)    -X ,   -Z ,    Y \t\t",
        " (21)     X ,   -Z ,   -Y \t\t(22)    -Y ,    X ,   -Z \t\t",
        " (23)    -Z ,   -Y ,    X \t\t(24)     Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "p 4/m m m": [
        " Space Group: P 4/m m m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "f m -3 c": [
        " Space Group: F m -3 c",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 192",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " (13)     Y ,    X ,1/2+Z \t(14) 1/2+Z ,    Y ,    X \t",
        " (15)     X ,1/2+Z ,    Y \t(16)     Y ,    X ,1/2-Z \t",
        " (17) 1/2-Z ,    Y ,    X \t(18)     X ,1/2-Z ,    Y \t",
        " (19) 1/2-Z ,    Y ,   -X \t(20)    -X ,1/2-Z ,    Y \t",
        " (21)     X ,1/2-Z ,   -Y \t(22)    -Y ,    X ,1/2-Z \t",
        " (23) 1/2-Z ,   -Y ,    X \t(24)     Y ,   -X ,1/2-Z \t",
        " ",
    ],
    "p n -3": [
        " Space Group: P n -3",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,1/2+Y ,   -Z \t",
        " ( 5)    -Z ,1/2+X ,1/2+Y \t( 6) 1/2+Y ,   -Z ,1/2+X \t",
        " ( 7) 1/2-Z ,    X ,1/2-Y \t( 8) 1/2-Y ,1/2-Z ,    X \t",
        " ( 9)     Y ,1/2-Z ,1/2-X \t(10) 1/2-X ,    Y ,1/2-Z \t",
        " (11) 1/2-Z ,1/2-X ,    Y \t(12)     X ,1/2-Y ,1/2-Z \t",
        " ",
    ],
    "p c c 2": [
        " Space Group: P c c 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,1/2+Z \t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "i 41 3 2": [
        " Space Group: I 41 3 2",
        " The lattice is noncentrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/4-Y ,3/4+X ,1/4+Z \t",
        " ( 5) 1/4+Z ,1/4-Y ,3/4+X \t( 6) 3/4+X ,1/4+Z ,1/4-Y \t",
        " ( 7) 1/2-X ,   -Y ,1/2+Z \t( 8)    -Z ,1/2+X ,1/2-Y \t",
        " ( 9) 1/2-Y ,   -Z ,1/2+X \t(10) 1/2+X ,1/2-Y ,   -Z \t",
        " (11) 1/2+Z ,1/2-X ,   -Y \t(12)    -Y ,1/2+Z ,1/2-X \t",
        " (13) 1/4+Y ,1/4-X ,3/4+Z \t(14) 3/4+Z ,1/4+Y ,1/4-X \t",
        " (15) 1/4-X ,3/4+Z ,1/4+Y \t(16) 3/4-X ,3/4-Z ,3/4-Y \t",
        " (17) 3/4-Y ,3/4-X ,3/4-Z \t(18) 3/4-Z ,3/4-Y ,3/4-X \t",
        " (19) 1/2+Y ,1/2-Z ,   -X \t(20) 3/4+Y ,1/4+X ,1/4-Z \t",
        " (21) 1/4-Z ,3/4+Y ,1/4+X \t(22) 1/4+X ,1/4-Z ,3/4+Y \t",
        " (23)    -X ,1/2+Y ,1/2-Z \t(24) 1/2-Z ,   -X ,1/2+Y \t",
        " ",
    ],
    "p 42 m c": [
        " Space Group: P 42 m c",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p 4 c c": [
        " Space Group: P 4 c c",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "p m -3 m": [
        " Space Group: P m -3 m",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " (13)     Y ,    X ,    Z \t\t(14)     Z ,    Y ,    X \t\t",
        " (15)     X ,    Z ,    Y \t\t(16)     Y ,    X ,   -Z \t\t",
        " (17)    -Z ,    Y ,    X \t\t(18)     X ,   -Z ,    Y \t\t",
        " (19)    -Z ,    Y ,   -X \t\t(20)    -X ,   -Z ,    Y \t\t",
        " (21)     X ,   -Z ,   -Y \t\t(22)    -Y ,    X ,   -Z \t\t",
        " (23)    -Z ,   -Y ,    X \t\t(24)     Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "p 32 1 2": [
        " Space Group: P 32 1 2",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 31m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,2/3+Z \t",
        " ( 3)    Y-X,   -X ,1/3+Z \t( 4)     X ,   X-Y,   -Z \t\t",
        " ( 5)    Y-X,    Y ,2/3-Z \t( 6)    -Y ,   -X ,1/3-Z \t",
        " ",
    ],
    "p 32 1 1": [
        " Space Group: P 32 1 1",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 3",
        " The Laue symmetry is 3",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,2/3+Z \t",
        " ( 3)    Y-X,   -X ,1/3+Z \t",
    ],
    "r -3 m r": [
        " Space Group: R -3 m r",
        " The lattice is centrosymmetric primitive rhombohedral",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 3mR",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     Y ,    X ,    Z \t\t",
        " ( 5)     Z ,    Y ,    X \t\t( 6)     X ,    Z ,    Y \t\t",
        " ",
    ],
    "p 3 c 1": [
        " Space Group: P 3 c 1",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3m1",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,1/2+Z \t",
        " ( 5)    -Y ,   -X ,1/2+Z \t( 6)     X ,   X-Y,1/2+Z \t",
        " ",
    ],
    "p 2 2 21": [
        " Space Group: P 2 2 21",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ( 3)    -X ,    Y ,1/2-Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 63": [
        " Space Group: P 63",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/2+Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,1/2+Z \t",
        " ",
    ],
    "p m 3": [
        " Space Group: P m 3",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " ",
    ],
    "p 42/m": [
        " Space Group: P 42/m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ",
    ],
    "p m c 21": [
        " Space Group: P m c 21",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,    Z \t\t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 42/n": [
        " Space Group: P 42/n",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,1/2+X ,1/2+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,    Z \t( 4) 1/2+Y ,   -X ,1/2+Z \t",
        " ",
    ],
    "a m a 2": [
        " Space Group: A m a 2",
        " The lattice is noncentrosymmetric A-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,    Z \t",
        " ( 3) 1/2+X ,   -Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 6/m": [
        " Space Group: P 6/m",
        " The lattice is centrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ",
    ],
    "p -6 c 2": [
        " Space Group: P -6 c 2",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    Y-X,   -X ,1/2-Z \t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)     X ,    Y ,1/2-Z \t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)    -Y ,   X-Y,1/2-Z \t",
        " ( 7)    Y-X,    Y ,1/2+Z \t( 8)     X ,   X-Y,   -Z \t\t",
        " ( 9)    -Y ,   -X ,1/2+Z \t(10)    Y-X,    Y ,   -Z \t\t",
        " (11)     X ,   X-Y,1/2+Z \t(12)    -Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "i -4 c 2": [
        " Space Group: I -4 c 2",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)     Y ,    X ,1/2-Z \t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)    -Y ,   -X ,1/2-Z \t",
        " ",
    ],
    "F -1": [
        " Space Group: F -1",
        " The lattice is centrosymmetric F-centered triclinic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is -1",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t",
    ],
    "p 3 1 m": [
        " Space Group: P 3 1 m",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 31m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     Y ,    X ,    Z \t\t",
        " ( 5)    -X ,   Y-X,    Z \t\t( 6)    X-Y,   -Y ,    Z \t\t",
        " ",
    ],
    "c c c 2": [
        " Space Group: C c c 2",
        " The lattice is noncentrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,1/2+Z \t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "i m 3": [
        " Space Group: I m 3",
        " The lattice is centrosymmetric I-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,    Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,    Y \t\t( 6)     Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,    X ,   -Y \t\t( 8)    -Y ,   -Z ,    X \t\t",
        " ( 9)     Y ,   -Z ,   -X \t\t(10)    -X ,    Y ,   -Z \t\t",
        " (11)    -Z ,   -X ,    Y \t\t(12)     X ,   -Y ,   -Z \t\t",
        " ",
    ],
    "p -4 3 m": [
        " Space Group: P -4 3 m",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     Y ,   -X ,   -Z \t\t",
        " ( 5)    -Z ,    Y ,   -X \t\t( 6)    -X ,   -Z ,    Y \t\t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13)    -Y ,    X ,   -Z \t\t(14)    -Z ,   -Y ,    X \t\t",
        " (15)     X ,   -Z ,   -Y \t\t(16)     X ,    Z ,    Y \t\t",
        " (17)     Y ,    X ,    Z \t\t(18)     Z ,    Y ,    X \t\t",
        " (19)     Y ,   -Z ,   -X \t\t(20)    -Y ,   -X ,    Z \t\t",
        " (21)     Z ,   -Y ,   -X \t\t(22)    -X ,    Z ,   -Y \t\t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "p -4 3 n": [
        " Space Group: P -4 3 n",
        " The lattice is noncentrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+Y ,1/2-X ,1/2-Z \t",
        " ( 5) 1/2-Z ,1/2+Y ,1/2-X \t( 6) 1/2-X ,1/2-Z ,1/2+Y \t",
        " ( 7)    -X ,   -Y ,    Z \t\t( 8)    -Z ,    X ,   -Y \t\t",
        " ( 9)    -Y ,   -Z ,    X \t\t(10)     X ,   -Y ,   -Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -Y ,    Z ,   -X \t\t",
        " (13) 1/2-Y ,1/2+X ,1/2-Z \t(14) 1/2-Z ,1/2-Y ,1/2+X \t",
        " (15) 1/2+X ,1/2-Z ,1/2-Y \t(16) 1/2+X ,1/2+Z ,1/2+Y \t",
        " (17) 1/2+Y ,1/2+X ,1/2+Z \t(18) 1/2+Z ,1/2+Y ,1/2+X \t",
        " (19)     Y ,   -Z ,   -X \t\t(20) 1/2-Y ,1/2-X ,1/2+Z \t",
        " (21) 1/2+Z ,1/2-Y ,1/2-X \t(22) 1/2-X ,1/2+Z ,1/2-Y \t",
        " (23)    -X ,    Y ,   -Z \t\t(24)    -Z ,   -X ,    Y \t\t",
        " ",
    ],
    "p -3 1 c": [
        " Space Group: P -3 1 c",
        " The lattice is centrosymmetric primitive trigonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 31m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     Y ,    X ,1/2+Z \t",
        " ( 5)    -X ,   Y-X,1/2+Z \t( 6)    X-Y,   -Y ,1/2+Z \t",
        " ",
    ],
    "r 3 m": [
        " Space Group: R 3 m",
        " The lattice is noncentrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 18",
        " The Laue symmetry is 3m1",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,    Z \t\t",
        " ( 5)    -Y ,   -X ,    Z \t\t( 6)     X ,   X-Y,    Z \t\t",
        " ",
    ],
    "p 21": [
        " Space Group: P 21",
        " The lattice is noncentrosymmetric primitive monoclinic",
        " Multiplicity of a general site is 2",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in y",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,   -Z \t",
        " ",
    ],
    "r -3": [
        " Space Group: R -3",
        " The lattice is centrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 18",
        " The Laue symmetry is 3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t",
    ],
    "c m": [
        " Space Group: C m",
        " The lattice is noncentrosymmetric C-centered monoclinic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is 2/m",
        " The unique monoclinic axis is b",
        " The location of the origin is arbitrary in x z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 32 2 1": [
        " Space Group: P 32 2 1",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3m1",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,2/3+Z \t",
        " ( 3)    Y-X,   -X ,1/3+Z \t( 4)     Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,   Y-X,2/3-Z \t( 6)    X-Y,   -Y ,1/3-Z \t",
        " ",
    ],
    "i 21 21 21": [
        " Space Group: I 21 21 21",
        " The lattice is noncentrosymmetric I-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2+X ,1/2-Y ,   -Z \t",
        " ( 3)    -X ,1/2+Y ,1/2-Z \t( 4) 1/2-X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 42 2 2": [
        " Space Group: P 42 2 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ( 5)    -X ,    Y ,   -Z \t\t( 6)    -Y ,   -X ,1/2-Z \t",
        " ( 7)     X ,   -Y ,   -Z \t\t( 8)     Y ,    X ,1/2-Z \t",
        " ",
    ],
    "i -4 2 m": [
        " Space Group: I -4 2 m",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,    Y ,   -Z \t\t( 6)     Y ,    X ,    Z \t\t",
        " ( 7)     X ,   -Y ,   -Z \t\t( 8)    -Y ,   -X ,    Z \t\t",
        " ",
    ],
    "p 65 1 1": [
        " Space Group: P 65 1 1",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,5/6+Z \t",
        " ( 3)    -Y ,   X-Y,2/3+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,1/3+Z \t( 6)     Y ,   Y-X,1/6+Z \t",
        " ",
    ],
    "p 61": [
        " Space Group: P 61",
        " The lattice is noncentrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 6/m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,1/6+Z \t",
        " ( 3)    -Y ,   X-Y,1/3+Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ( 5)    Y-X,   -X ,2/3+Z \t( 6)     Y ,   Y-X,5/6+Z \t",
        " ",
    ],
    "i 2 3": [
        " Space Group: I 2 3",
        " The lattice is noncentrosymmetric I-centered cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,   -Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,   -Y \t\t( 6)    -Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,   -X ,    Y \t\t( 8)     Y ,   -Z ,   -X \t\t",
        " ( 9)    -Y ,    Z ,   -X \t\t(10)    -X ,   -Y ,    Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "i -4 2 d": [
        " Space Group: I -4 2 d",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5) 1/2-X ,    Y ,3/4-Z \t( 6)     Y ,1/2+X ,1/4+Z \t",
        " ( 7) 1/2+X ,   -Y ,3/4-Z \t( 8)    -Y ,1/2-X ,1/4+Z \t",
        " ",
    ],
    "p a 3": [
        " Space Group: P a 3",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 24",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,    Y ,1/2-Z \t",
        " ( 5) 1/2-Z ,1/2+X ,    Y \t( 6)     Y ,1/2-Z ,1/2+X \t",
        " ( 7)    -Z ,1/2+X ,1/2-Y \t( 8) 1/2-Y ,   -Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/2-Z ,   -X \t(10)    -X ,1/2+Y ,1/2-Z \t",
        " (11) 1/2-Z ,   -X ,1/2+Y \t(12) 1/2+X ,1/2-Y ,   -Z \t",
        " ",
    ],
    "f 2 3": [
        " Space Group: F 2 3",
        " The lattice is noncentrosymmetric F-centered cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     X ,   -Y ,   -Z \t\t",
        " ( 5)    -Z ,    X ,   -Y \t\t( 6)    -Y ,   -Z ,    X \t\t",
        " ( 7)    -Z ,   -X ,    Y \t\t( 8)     Y ,   -Z ,   -X \t\t",
        " ( 9)    -Y ,    Z ,   -X \t\t(10)    -X ,   -Y ,    Z \t\t",
        " (11)     Z ,   -X ,   -Y \t\t(12)    -X ,    Y ,   -Z \t\t",
        " ",
    ],
    "i 4 c m": [
        " Space Group: I 4 c m",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)    -Y ,   -X ,1/2+Z \t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)     Y ,    X ,1/2+Z \t",
        " ",
    ],
    "r 3 c": [
        " Space Group: R 3 c",
        " The lattice is noncentrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 18",
        " The Laue symmetry is 3m1",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,1/2+Z \t",
        " ( 5)    -Y ,   -X ,1/2+Z \t( 6)     X ,   X-Y,1/2+Z \t",
        " ",
    ],
    "p n m a": [
        " Space Group: P n m a",
        " The lattice is centrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,1/2+Y ,1/2+Z \t",
        " ( 3)     X ,1/2-Y ,    Z \t( 4) 1/2-X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "r 3 c r": [
        " Space Group: R 3 c r",
        " The lattice is noncentrosymmetric primitive rhombohedral",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3mR",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+Y ,1/2+X ,1/2+Z \t",
        " ( 5) 1/2+Z ,1/2+Y ,1/2+X \t( 6) 1/2+X ,1/2+Z ,1/2+Y \t",
        " ",
    ],
    "p n c 2": [
        " Space Group: P n c 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,1/2+Z \t",
        " ( 3)     X ,1/2-Y ,1/2+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "c 2 2 21": [
        " Space Group: C 2 2 21",
        " The lattice is noncentrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ( 3)    -X ,    Y ,1/2-Z \t( 4)    -X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "r 3 m r": [
        " Space Group: R 3 m r",
        " The lattice is noncentrosymmetric primitive rhombohedral",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 3mR",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4)     Y ,    X ,    Z \t\t",
        " ( 5)     Z ,    Y ,    X \t\t( 6)     X ,    Z ,    Y \t\t",
        " ",
    ],
    "p 43 2 2": [
        " Space Group: P 43 2 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,3/4+Z \t",
        " ( 3)    -X ,   -Y ,1/2+Z \t( 4)     Y ,   -X ,1/4+Z \t",
        " ( 5)    -X ,    Y ,   -Z \t\t( 6)    -Y ,   -X ,3/4-Z \t",
        " ( 7)     X ,   -Y ,1/2-Z \t( 8)     Y ,    X ,1/4-Z \t",
        " ",
    ],
    "r 3 2": [
        " Space Group: R 3 2",
        " The lattice is noncentrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 18",
        " The Laue symmetry is 3m1",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     Y ,    X ,   -Z \t\t",
        " ( 5)    -X ,   Y-X,   -Z \t\t( 6)    X-Y,   -Y ,   -Z \t\t",
        " ",
    ],
    "p m a 2": [
        " Space Group: P m a 2",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,    Z \t",
        " ( 3) 1/2+X ,   -Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "i 4/m m m": [
        " Space Group: I 4/m m m",
        " The lattice is centrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 32",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,    Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,    Z \t\t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7)     X ,   -Y ,    Z \t\t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "c c c a": [
        " Space Group: C c c a",
        " The lattice is centrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-X ,    Y ,1/2+Z \t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4) 1/2-X ,   -Y ,    Z \t",
        " ",
    ],
    "i 41 m d": [
        " Space Group: I 41 m d",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,1/2+X ,1/4+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,1/2+Z \t( 4) 1/2+Y ,   -X ,3/4+Z \t",
        " ( 5)    -X ,    Y ,    Z \t\t( 6)    -Y ,1/2-X ,1/4+Z \t",
        " ( 7) 1/2+X ,1/2-Y ,1/2+Z \t( 8) 1/2+Y ,    X ,3/4+Z \t",
        " ",
    ],
    "c c c m": [
        " Space Group: C c c m",
        " The lattice is centrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,    Y ,1/2+Z \t",
        " ( 3)     X ,   -Y ,1/2+Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p 41 21 2": [
        " Space Group: P 41 21 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,1/2+X ,1/4+Z \t",
        " ( 3)    -X ,   -Y ,1/2+Z \t( 4) 1/2+Y ,1/2-X ,3/4+Z \t",
        " ( 5) 1/2-X ,1/2+Y ,1/4-Z \t( 6)    -Y ,   -X ,1/2-Z \t",
        " ( 7) 1/2+X ,1/2-Y ,3/4-Z \t( 8)     Y ,    X ,   -Z \t\t",
        " ",
    ],
    "p 31": [
        " Space Group: P 31",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 3",
        " The Laue symmetry is 3",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,1/3+Z \t",
        " ( 3)    Y-X,   -X ,2/3+Z \t",
    ],
    "p 32": [
        " Space Group: P 32",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 3",
        " The Laue symmetry is 3",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,2/3+Z \t",
        " ( 3)    Y-X,   -X ,1/3+Z \t",
    ],
    "p 42/m n m": [
        " Space Group: P 42/m n m",
        " The lattice is centrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2-Y ,1/2+X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4) 1/2+Y ,1/2-X ,1/2+Z \t",
        " ( 5) 1/2-X ,1/2+Y ,1/2+Z \t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7) 1/2+X ,1/2-Y ,1/2+Z \t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p 3 1 2": [
        " Space Group: P 3 1 2",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 31m",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     X ,   X-Y,   -Z \t\t",
        " ( 5)    Y-X,    Y ,   -Z \t\t( 6)    -Y ,   -X ,   -Z \t\t",
        " ",
    ],
    "i 41 2 2": [
        " Space Group: I 41 2 2",
        " The lattice is noncentrosymmetric I-centered tetragonal",
        " Multiplicity of a general site is 16",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,1/2+X ,1/4+Z \t",
        " ( 3) 1/2-X ,1/2-Y ,1/2+Z \t( 4) 1/2+Y ,   -X ,3/4+Z \t",
        " ( 5) 1/2-X ,    Y ,3/4-Z \t( 6)    -Y ,   -X ,   -Z \t\t",
        " ( 7)     X ,1/2-Y ,1/4-Z \t( 8) 1/2+Y ,1/2+X ,1/2-Z \t",
        " ",
    ],
    "p -3 m 1": [
        " Space Group: P -3 m 1",
        " The lattice is centrosymmetric primitive trigonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 3m1",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)    Y-X,    Y ,    Z \t\t",
        " ( 5)    -Y ,   -X ,    Z \t\t( 6)     X ,   X-Y,    Z \t\t",
        " ",
    ],
    "a b m 2": [
        " Space Group: A b m 2",
        " The lattice is noncentrosymmetric A-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -X ,1/2+Y ,    Z \t",
        " ( 3)     X ,1/2-Y ,    Z \t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p n -3 n": [
        " Space Group: P n -3 n",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,1/2+Y ,   -Z \t",
        " ( 5)    -Z ,1/2+X ,1/2+Y \t( 6) 1/2+Y ,   -Z ,1/2+X \t",
        " ( 7) 1/2-Z ,    X ,1/2-Y \t( 8) 1/2-Y ,1/2-Z ,    X \t",
        " ( 9)     Y ,1/2-Z ,1/2-X \t(10) 1/2-X ,    Y ,1/2-Z \t",
        " (11) 1/2-Z ,1/2-X ,    Y \t(12)     X ,1/2-Y ,1/2-Z \t",
        " (13) 1/2+Y ,1/2+X ,1/2+Z \t(14) 1/2+Z ,1/2+Y ,1/2+X \t",
        " (15) 1/2+X ,1/2+Z ,1/2+Y \t(16)     Y ,    X ,1/2-Z \t",
        " (17) 1/2-Z ,    Y ,    X \t(18)     X ,1/2-Z ,    Y \t",
        " (19)    -Z ,1/2+Y ,   -X \t(20)    -X ,   -Z ,1/2+Y \t",
        " (21) 1/2+X ,   -Z ,   -Y \t(22)    -Y ,1/2+X ,   -Z \t",
        " (23)    -Z ,   -Y ,1/2+X \t(24) 1/2+Y ,   -X ,   -Z \t",
        " ",
    ],
    "r 3": [
        " Space Group: R 3",
        " The lattice is noncentrosymmetric R-centered trigonal",
        " Multiplicity of a general site is 9",
        " The Laue symmetry is 3",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t",
    ],
    "c 2 2 2": [
        " Space Group: C 2 2 2",
        " The lattice is noncentrosymmetric C-centered orthorhombic",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        "\n (0,0,0; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     X ,   -Y ,   -Z \t\t",
        " ( 3)    -X ,    Y ,   -Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ",
    ],
    "p n -3 m": [
        " Space Group: P n -3 m",
        " The lattice is centrosymmetric primitive cubic",
        " Multiplicity of a general site is 48",
        " The Laue symmetry is m3m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/2+X ,1/2+Y ,   -Z \t",
        " ( 5)    -Z ,1/2+X ,1/2+Y \t( 6) 1/2+Y ,   -Z ,1/2+X \t",
        " ( 7) 1/2-Z ,    X ,1/2-Y \t( 8) 1/2-Y ,1/2-Z ,    X \t",
        " ( 9)     Y ,1/2-Z ,1/2-X \t(10) 1/2-X ,    Y ,1/2-Z \t",
        " (11) 1/2-Z ,1/2-X ,    Y \t(12)     X ,1/2-Y ,1/2-Z \t",
        " (13)     Y ,    X ,    Z \t\t(14)     Z ,    Y ,    X \t\t",
        " (15)     X ,    Z ,    Y \t\t(16) 1/2+Y ,1/2+X ,   -Z \t",
        " (17)    -Z ,1/2+Y ,1/2+X \t(18) 1/2+X ,   -Z ,1/2+Y \t",
        " (19) 1/2-Z ,    Y ,1/2-X \t(20) 1/2-X ,1/2-Z ,    Y \t",
        " (21)     X ,1/2-Z ,1/2-Y \t(22) 1/2-Y ,    X ,1/2-Z \t",
        " (23) 1/2-Z ,1/2-Y ,    X \t(24)     Y ,1/2-X ,1/2-Z \t",
        " ",
    ],
    "p 42 c m": [
        " Space Group: P 42 c m",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/2+Z \t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)     Y ,   -X ,1/2+Z \t",
        " ( 5)    -X ,    Y ,1/2+Z \t( 6)    -Y ,   -X ,    Z \t\t",
        " ( 7)     X ,   -Y ,1/2+Z \t( 8)     Y ,    X ,    Z \t\t",
        " ",
    ],
    "p 6/m 1 1": [
        " Space Group: P 6/m 1 1",
        " The lattice is centrosymmetric primitive hexagonal",
        " Multiplicity of a general site is 12",
        " The Laue symmetry is 6/m",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    X-Y,    X ,    Z \t\t",
        " ( 3)    -Y ,   X-Y,    Z \t\t( 4)    -X ,   -Y ,    Z \t\t",
        " ( 5)    Y-X,   -X ,    Z \t\t( 6)     Y ,   Y-X,    Z \t\t",
        " ",
    ],
    "p 21 21 21": [
        " Space Group: P 21 21 21",
        " The lattice is noncentrosymmetric primitive orthorhombic",
        " Multiplicity of a general site is 4",
        " The Laue symmetry is mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2) 1/2+X ,1/2-Y ,   -Z \t",
        " ( 3)    -X ,1/2+Y ,1/2-Z \t( 4) 1/2-X ,   -Y ,1/2+Z \t",
        " ",
    ],
    "f d -3": [
        " Space Group: F d -3",
        " The lattice is centrosymmetric F-centered cubic",
        " Multiplicity of a general site is 96",
        " The Laue symmetry is m3",
        " The inversion center is located at 0,0,0",
        "\n The equivalent positions are:",
        "\n (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+\n",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Z ,    X ,    Y \t\t",
        " ( 3)     Y ,    Z ,    X \t\t( 4) 1/4+X ,1/4+Y ,   -Z \t",
        " ( 5)    -Z ,1/4+X ,1/4+Y \t( 6) 1/4+Y ,   -Z ,1/4+X \t",
        " ( 7) 1/4-Z ,1/2+X ,3/4-Y \t( 8) 3/4-Y ,1/4-Z ,1/2+X \t",
        " ( 9) 1/2+Y ,1/4-Z ,3/4-X \t(10) 3/4-X ,1/2+Y ,1/4-Z \t",
        " (11) 1/4-Z ,3/4-X ,1/2+Y \t(12) 1/2+X ,3/4-Y ,1/4-Z \t",
        " ",
    ],
    "p -4 b 2": [
        " Space Group: P -4 b 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)     Y ,   -X ,   -Z \t\t",
        " ( 3)    -X ,   -Y ,    Z \t\t( 4)    -Y ,    X ,   -Z \t\t",
        " ( 5) 1/2-X ,1/2+Y ,    Z \t( 6) 1/2+Y ,1/2+X ,   -Z \t",
        " ( 7) 1/2+X ,1/2-Y ,    Z \t( 8) 1/2-Y ,1/2-X ,   -Z \t",
        " ",
    ],
    "p 3 1 c": [
        " Space Group: P 3 1 c",
        " The lattice is noncentrosymmetric primitive trigonal",
        " Multiplicity of a general site is 6",
        " The Laue symmetry is 31m",
        " The location of the origin is arbitrary in z",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,   X-Y,    Z \t\t",
        " ( 3)    Y-X,   -X ,    Z \t\t( 4)     Y ,    X ,1/2+Z \t",
        " ( 5)    -X ,   Y-X,1/2+Z \t( 6)    X-Y,   -Y ,1/2+Z \t",
        " ",
    ],
    "p 41 2 2": [
        " Space Group: P 41 2 2",
        " The lattice is noncentrosymmetric primitive tetragonal",
        " Multiplicity of a general site is 8",
        " The Laue symmetry is 4/mmm",
        "\n The equivalent positions are:",
        " ( 1)     X ,    Y ,    Z \t\t( 2)    -Y ,    X ,1/4+Z \t",
        " ( 3)    -X ,   -Y ,1/2+Z \t( 4)     Y ,   -X ,3/4+Z \t",
        " ( 5)    -X ,    Y ,   -Z \t\t( 6)    -Y ,   -X ,1/4-Z \t",
        " ( 7)     X ,   -Y ,1/2-Z \t( 8)     Y ,    X ,3/4-Z \t",
        " ",
    ],
}
