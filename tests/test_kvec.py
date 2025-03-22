import os
import sys
import numpy as np
import pytest

import importlib.util
if importlib.util.find_spec('GSASII') is None:
    home = os.path.dirname(__file__)
    sys.path.append(os.path.dirname(home))

import GSASII.k_vector_search as kvs

gen_pos = False

pcell = [
    [5.144191, -2.970000, -0.000000],
    [0.000000, 5.940000, -0.000000],
    [0.000000, 0.000000, 16.700001]
]

ppos = [
    [0.000000, 0.000000, 0.000000],
    [0.000000, 0.000000, 0.500000],
    [0.000000, 0.000000, 0.330000],
    [0.000000, 0.000000, 0.670000],
    [0.333333, 0.666667, 0.666667],
    [0.333333, 0.666667, 0.166667],
    [0.333333, 0.666667, 0.996667],
    [0.333333, 0.666667, 0.336667],
    [0.666667, 0.333333, 0.333333],
    [0.666667, 0.333333, 0.833333],
    [0.666667, 0.333333, 0.663333],
    [0.666667, 0.333333, 0.003333],
    [0.333000, 0.000000, 0.251000],
    [0.000000, 0.333000, 0.251000],
    [0.667000, 0.667000, 0.251000],
    [0.667000, 0.000000, 0.749000],
    [0.000000, 0.667000, 0.749000],
    [0.333000, 0.333000, 0.749000],
    [0.666333, 0.666667, 0.917667],
    [0.333333, 0.999667, 0.917667],
    [0.000333, 0.333667, 0.917667],
    [0.000333, 0.666667, 0.415667],
    [0.333333, 0.333667, 0.415667],
    [0.666333, 0.999667, 0.415667],
    [0.999667, 0.333333, 0.584333],
    [0.666667, 0.666333, 0.584333],
    [0.333667, 0.000333, 0.584333],
    [0.333667, 0.333333, 0.082333],
    [0.666667, 0.000333, 0.082333],
    [0.999667, 0.666333, 0.082333]
]

cr_atoms = [24 for _ in range(12)]
s_atoms = [16 for _ in range(18)]
nums = cr_atoms + s_atoms

nuc_p = list()
for i in range(6):
    for j in range(6):
        for k in range(6):
            nuc_p.append([i, j, k])

spos = [
    5.566667,
    4.916235,
    4.379751
]

threshold = 0.008
brav_type = "hR"

try:
    import seekpath
    k_search = kvs.kVector(
        brav_type,
        pcell,
        ppos,
        nums,
        nuc_p,
        spos,
        threshold
    )
except ModuleNotFoundError:
    k_search = None

selftestlist = []
selftestquiet = True


def _ReportTest():
    '''Report name and doc string of current routine when ``selftestquiet``
    is False'
    '''
    if not selftestquiet:
        import inspect
        caller = inspect.stack()[1][3]
        doc = eval(caller).__doc__
        if doc is not None:
            print(
                f'testing {os.path.split(__file__)[1]}'
                f' with {caller}\n\t({doc})'
            )
        else:
            print(f'testing {os.path.split(__file__)[1]} with {caller}')


def test_AtomID():
    '''self-test #0: test the unique ID generation routine'''
    _ReportTest()
    atom_types = ["Cr", "Cr", "Cr", "Sb", "Sb", "Se"]
    atom_types_id = kvs.unique_id_gen(atom_types)
    msg = "Unique ID generation failed"
    assert atom_types_id == [1, 1, 1, 2, 2, 3], msg
    print("test_AtomID passed")


selftestlist.append(test_AtomID)


@pytest.mark.skipif(k_search is None, reason='No seekpath module')
def test_LatConstruct():
    '''self-test #1: test the lattice vectors construction routine'''
    _ReportTest()
    cell_params = [5., 5., 5., 90, 90., 90.]
    latt_vec = kvs.lat_params_to_vec(cell_params)
    msg = "Lattice vectors construction failed"
    assert np.allclose(latt_vec[0], [5., 0., 0.]), msg
    assert np.allclose(latt_vec[1], [0., 5., 0.]), msg
    assert np.allclose(latt_vec[2], [0., 0., 5.]), msg

    cell_params = [
        5., 5., 10.,
        90., 90., 120.
    ]
    latt_vec = kvs.lat_params_to_vec(cell_params)
    assert np.allclose(
        latt_vec[0], [4.33013, -2.5, 0.],
        atol=1e-5
    ), msg
    assert np.allclose(
        latt_vec[1], [0., 5.0, 0.],
        atol=1e-5
    ), msg
    assert np.allclose(
        latt_vec[2], [0., 0., 10.],
        atol=1e-5
    ), msg

    cell_params = [
        5., 6., 7.,
        91., 95., 112.
    ]
    latt_vec = kvs.lat_params_to_vec(cell_params)
    assert np.allclose(
        latt_vec[0], [4.61218, -1.88092, -0.43578],
        atol=1e-5
    ), msg
    assert np.allclose(
        latt_vec[1], [0., 5.99909, -0.10471],
        atol=1e-5
    ), msg
    assert np.allclose(
        latt_vec[2], [0., 0., 7.],
        atol=1e-5
    ), msg

    print("test_LatConstruct passed")


selftestlist.append(test_LatConstruct)


@pytest.mark.skipif(k_search is None, reason='No seekpath module')
def test_CriticalRoutines():
    '''self-test #2: test the critical routines'''
    _ReportTest()
    hkl_p = k_search.hklConvToPrim(nuc_p[2][:3])
    msg = "Test for conventional to primitive setting failed"
    assert np.allclose(
        hkl_p,
        [0.66667, 0.66667, 0.66667],
        atol=1e-5
    ), msg
    msg = "Test for primitive to conventional setting failed"
    k_conv = k_search.kVecPrimToConv([0.5, 0.5, 0])
    assert np.allclose(
        k_conv,
        [0., .5, 1.],
        atol=1e-5
    ), msg


selftestlist.append(test_CriticalRoutines)


@pytest.mark.skipif(k_search is None, reason='No seekpath module')
def test_KVecCandidateUpdate():
    '''self-test #3: test the updating of the list of alternative k vectors'''
    _ReportTest()

    msg = "Test for updating the list of alternative k vectors failed"

    k_trial = [.5, .5, .5]  # T point
    k_opt_list = list()
    k_opt_dist = list()
    k_opt_ad = list()
    k_opt_md = list()
    k_opt_out = k_search.updateCandidateList(
        k_trial,
        k_opt_list,
        k_opt_dist,
        k_opt_ad,
        k_opt_md,
        False
    )
    (k_opt_list, k_opt_dist, k_opt_ad, k_opt_md) = k_opt_out

    k_trial = [.5, 0, .5]  # F point
    k_opt_out = k_search.updateCandidateList(
        k_trial,
        k_opt_list,
        k_opt_dist,
        k_opt_ad,
        k_opt_md,
        False
    )
    (k_opt_list, k_opt_dist, k_opt_ad, k_opt_md) = k_opt_out

    assert np.allclose(
        k_opt_list[0],
        [.5, 0., .5],
        atol=1e-5
    ), msg
    assert np.allclose(
        k_opt_dist[0],
        0.14445,
        atol=1e-5
    ), msg

    assert np.allclose(
        k_opt_list[1],
        [.5, .5, .5],
        atol=1e-5
    ), msg
    assert np.allclose(
        k_opt_dist[1],
        .20353,
        atol=1e-5
    ), msg

    k_trial = [0, 0, 0]  # Gamma point
    k_opt_out = k_search.updateCandidateList(
        k_trial,
        k_opt_list,
        k_opt_dist,
        k_opt_ad,
        k_opt_md,
        False
    )
    (k_opt_list, k_opt_dist, k_opt_ad, k_opt_md) = k_opt_out

    assert np.allclose(
        k_opt_list[0],
        [0., 0., 0.],
        atol=1e-5
    ), msg
    assert np.allclose(
        k_opt_dist[0],
        0.,
        atol=1e-5
    ), msg

    k_opt_final, kd_opt_final, _, _ = k_search.kOptFinder()
    assert np.allclose(
        k_opt_final[0],
        [0., 0., 0.],
        atol=1e-5
    ), msg
    assert np.allclose(
        kd_opt_final[0],
        0.,
        atol=1e-5
    ), msg


selftestlist.append(test_KVecCandidateUpdate)


@pytest.mark.skipif(k_search is None, reason='No seekpath module')
def test_KVecSearch():
    '''self-test #4: test the k vector search routine'''
    _ReportTest()

    # High symmetry point.
    lambda_val = 2.4109
    two_thetas = [
        24.613175764104273,
        35.08694489435111,
        56.92736582327283,
        62.945222217048965,
        79.49838666374927,
        84.75564148306853,
        100.43822871252318,
        105.78539737883935,
        122.99828511273124,
        129.45511757798127
    ]

    spos_gen = list()
    for two_theta in two_thetas:
        d_tmp = lambda_val / (2. * np.sin(np.deg2rad(two_theta) / 2.))
        spos_gen.append(d_tmp)

    latt_params = [
        5.655600,
        5.655600,
        5.655600,
        90.000000,
        90.000000,
        90.000000
    ]
    latt_vec = kvs.lat_params_to_vec(latt_params)

    atoms_labels = [1, 1, 1, 1]

    atom_pos_p1 = [
        [0.000000, 0.000000, 0.000000],
        [0.000000, 0.500000, 0.500000],
        [0.500000, 0.000000, 0.500000],
        [0.500000, 0.500000, 0.000000]
    ]

    hkl_refls = list()
    for i in range(6):
        for j in range(6):
            for k in range(6):
                hkl_refls.append([i, j, k])

    brav_type = "cF"
    threshold = 1.E-6

    k_search = kvs.kVector(
        brav_type,
        latt_vec,
        atom_pos_p1,
        atoms_labels,
        hkl_refls,
        spos_gen,
        threshold,
        option=2,
        kstep=[0.002, 0.002, 0.002],
        processes=16
    )

    k_opt = k_search.kOptFinder()
    k_opt_final = k_search.kVecPrimToConv(k_opt[0])

    msg = "Test for k vector search failed"
    assert np.allclose(
        k_opt_final[0],
        [0., 1., 0.],
        atol=1e-5
    ), msg
    assert np.allclose(
        k_opt[1][0],
        0.,
        atol=1e-5
    ), msg

    # High symmetry direction.
    lambda_val = 2.4109
    two_thetas = [
        3.989827234732058,
        24.694684773857567,
        29.27168183349248,
        31.061600212829333,
        31.24029252242816,
        34.877431565785194,
        40.21710077021546,
        46.02390422155341,
        54.92232759586314,
        62.073039736321064
    ]

    spos_gen = list()
    for two_theta in two_thetas:
        d_tmp = lambda_val / (2. * np.sin(np.deg2rad(two_theta) / 2.))
        spos_gen.append(d_tmp)

    latt_params = [
        5.598000,
        5.598000,
        8.955600,
        90.000000,
        90.000000,
        120.000000
    ]
    latt_vec = kvs.lat_params_to_vec(latt_params)

    atoms_labels = [1]

    atom_pos_p1 = [
        [0.000000, 0.000000, 0.000000]
    ]

    hkl_refls = list()
    for i in range(6):
        for j in range(6):
            for k in range(6):
                hkl_refls.append([i, j, k])

    brav_type = "hP"
    threshold = 1.E-6

    k_search = kvs.kVector(
        brav_type,
        latt_vec,
        atom_pos_p1,
        atoms_labels,
        hkl_refls,
        spos_gen,
        threshold,
        option=2,
        kstep=[0.002, 0.002, 0.002],
        processes=16
    )

    k_opt = k_search.kOptFinder()
    k_opt_final = k_search.kVecPrimToConv(k_opt[0])

    assert np.allclose(
        k_opt_final[0],
        [.14, 0., 0.],
        atol=1e-5
    ), msg
    assert np.allclose(
        k_opt[1][0],
        0.,
        atol=1e-5
    ), msg


selftestlist.append(test_KVecSearch)


if __name__ == "__main__":
    selftestquiet = False
    for test in selftestlist:
        test()
    print("All k-vector tests passed.")
