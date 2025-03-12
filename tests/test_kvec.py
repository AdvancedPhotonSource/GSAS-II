import GSASII.k_vector_search as kvs
import numpy as np
import json

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

# each type of atoms should have its own unique
# integer label. This is necessary for the symmetry
# identification from the input structure.
cr_atoms = [24 for _ in range(12)]
s_atoms = [16 for _ in range(18)]
nums = cr_atoms + s_atoms

# nucleus peaks
nuc_p = list()
for i in range(6):
    for j in range(6):
        for k in range(6):
            nuc_p.append([i, j, k])

# satellite peaks, in d-spacing
spos = [
    5.566667,
    4.916235,
    4.379751
]

# criterion for determining whether a unique k vector
# can be determined. The delta_d/d characteristic value
# for the insturment resolution should be used here.
threshold = 0.008

# the Bravais lattice type of the nucleus structure
brav_type = "hR"
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#               ^^^^^ Inputs ^^^^^
#               |||||        |||||
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

if __name__ == "__main__":
    k_search = kvs.kVector(
        brav_type,
        pcell,
        ppos,
        nums,
        nuc_p,
        spos,
        threshold
    )

    # test the unique ID generation routine
    atom_types = ["Cr", "Cr", "Cr", "Sb", "Sb", "Se"]
    atom_types_id = kvs.unique_id_gen(atom_types)
    print("Unique ID generation\n===")
    print("Atom types present: ", atom_types)
    print("Atom type IDs: ", atom_types_id)

    # test the lattice vectors construction routine
    #
    # cubic lattice
    cell_params = [5., 5., 5., 90, 90., 90.]
    latt_vec = kvs.lat_params_to_vec(cell_params)
    print("\nLattice vectors construction\n===")
    print("Cell parameters: ", cell_params)
    print("Lattice vectors: ", latt_vec)
    print("|a| = ", np.linalg.norm(np.array(latt_vec[0])))
    print("|b| = ", np.linalg.norm(np.array(latt_vec[1])))
    print("|c| = ", np.linalg.norm(np.array(latt_vec[2])))
    bc_dot_P = np.dot(np.array(latt_vec[1]), np.array(latt_vec[2]))
    ac_dot_P = np.dot(np.array(latt_vec[0]), np.array(latt_vec[2]))
    ab_dot_P = np.dot(np.array(latt_vec[0]), np.array(latt_vec[1]))
    a_norm = np.linalg.norm(np.array(latt_vec[0]))
    b_norm = np.linalg.norm(np.array(latt_vec[1]))
    c_norm = np.linalg.norm(np.array(latt_vec[2]))
    print("alpha = ", np.rad2deg(np.arccos(bc_dot_P / (b_norm * c_norm))))
    print("beta = ", np.rad2deg(np.arccos(ac_dot_P / (a_norm * c_norm))))
    print("gamma = ", np.rad2deg(np.arccos(ab_dot_P / (a_norm * b_norm))))

    # rhombohedral lattice
    cell_params = [
        5.,
        5.,
        10.,
        90.,
        90.,
        120.
    ]
    latt_vec = kvs.lat_params_to_vec(cell_params)
    print("\nCell parameters: ", cell_params)
    print("Lattice vectors: ", latt_vec)
    print("|a| = ", np.linalg.norm(np.array(latt_vec[0])))
    print("|b| = ", np.linalg.norm(np.array(latt_vec[1])))
    print("|c| = ", np.linalg.norm(np.array(latt_vec[2])))
    bc_dot_P = np.dot(np.array(latt_vec[1]), np.array(latt_vec[2]))
    ac_dot_P = np.dot(np.array(latt_vec[0]), np.array(latt_vec[2]))
    ab_dot_P = np.dot(np.array(latt_vec[0]), np.array(latt_vec[1]))
    a_norm = np.linalg.norm(np.array(latt_vec[0]))
    b_norm = np.linalg.norm(np.array(latt_vec[1]))
    c_norm = np.linalg.norm(np.array(latt_vec[2]))
    print("alpha = ", np.rad2deg(np.arccos(bc_dot_P / (b_norm * c_norm))))
    print("beta = ", np.rad2deg(np.arccos(ac_dot_P / (a_norm * c_norm))))
    print("gamma = ", np.rad2deg(np.arccos(ab_dot_P / (a_norm * b_norm))))

    # triclinic lattice
    cell_params = [
        5.,
        6.,
        7.,
        91.,
        95.,
        112.
    ]
    latt_vec = kvs.lat_params_to_vec(cell_params)
    print("\nCell parameters: ", cell_params)
    print("Lattice vectors: ", latt_vec)
    print("|a| = ", np.linalg.norm(np.array(latt_vec[0])))
    print("|b| = ", np.linalg.norm(np.array(latt_vec[1])))
    print("|c| = ", np.linalg.norm(np.array(latt_vec[2])))
    bc_dot_P = np.dot(np.array(latt_vec[1]), np.array(latt_vec[2]))
    ac_dot_P = np.dot(np.array(latt_vec[0]), np.array(latt_vec[2]))
    ab_dot_P = np.dot(np.array(latt_vec[0]), np.array(latt_vec[1]))
    a_norm = np.linalg.norm(np.array(latt_vec[0]))
    b_norm = np.linalg.norm(np.array(latt_vec[1]))
    c_norm = np.linalg.norm(np.array(latt_vec[2]))
    print("alpha = ", np.rad2deg(np.arccos(bc_dot_P / (b_norm * c_norm))))
    print("beta = ", np.rad2deg(np.arccos(ac_dot_P / (a_norm * c_norm))))
    print("gamma = ", np.rad2deg(np.arccos(ab_dot_P / (a_norm * b_norm))))

    # test the k-path search engine
    k_path = k_search.kpathFinder()
    formatted_json = json.dumps(k_path, indent=4)
    print("\nk search path\n===")
    print(formatted_json)

    # test the conversion for hkl from the conventional to
    # the primitive setting
    hkl_p = k_search.hklConvToPrim(nuc_p[2][:3])
    print("\nPrimitive hkl for (101)\n===")
    print(hkl_p)

    # test the conversion of the k vector from the primitive
    # to the conventional setting
    k_conv = k_search.kVecPrimToConv([0.5, 0.5, 0])
    print("\nConventional k vector for (0.5, 0.5, 0)\n===")
    print(k_conv)

    # test the generation of a k point along a certain
    # vector
    k_point = k_search.pointOnVector(
        [0, 0, 0],
        [.5, .5, .5],
        0.433
    )
    print("\nk point in the middle of (0, 0, 0) -> (.5, .5, .5)\n===")
    print(k_point)

    # test the insersion of a value while keeping the
    # order of entries in the list
    test_list = [1, 2, 3, 5]
    val_tmp = 4
    new_list = k_search.insIntoSortedList(test_list, val_tmp)
    print("\nNew list:\n===")
    print(new_list[0])
    print(f"Position of '{val_tmp}' in new list\n===")
    print(new_list[1])

    # test the updating of the list of alternative k vectors,
    # given a trial k vector of (.5, .5, .5)
    k_trial = [.5, .5, .5]  # T point
    k_opt_list = list()
    k_opt_dist = list()
    k_opt_out = k_search.updateCandidateList(
        k_trial,
        k_opt_list,
        k_opt_dist,
        False
    )
    (k_opt_list, k_opt_dist) = k_opt_out
    print("\nOptimized k vector alternatives\n===")
    print(k_opt_list)
    print("Indicator distances of alternative k vectors\n===")
    print(k_opt_dist)

    # test the updating of the list of alternative k vectors,
    # given a trial k vector of (.5, 0, .5)
    k_trial = [.5, 0, .5]  # F point
    k_opt_out = k_search.updateCandidateList(
        k_trial,
        k_opt_list,
        k_opt_dist,
        False
    )
    (k_opt_list, k_opt_dist) = k_opt_out
    print("\nOptimized k vector alternatives\n===")
    print(k_opt_list)
    print("Indicator distances of alternative k vectors\n===")
    print(k_opt_dist)

    # test the updating of the list of alternative k vectors,
    # given the trial k vector of (0, 0, 0)
    k_trial = [0, 0, 0]  # Gamma point
    k_opt_out = k_search.updateCandidateList(
        k_trial,
        k_opt_list,
        k_opt_dist,
        False
    )
    (k_opt_list, k_opt_dist) = k_opt_out
    print("\nOptimized k vector alternatives\n===")
    print(k_opt_list)
    print("Indicator distances of alternative k vectors\n===")
    print(k_opt_dist)

    # test the overall k vector search routine, putting
    # together all the little pieces tested above
    k_opt_final = k_search.kOptFinder()
    print("\nOptimal candidate of k vector\n===")
    print(k_opt_final)

    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    # High symmetry point.
    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    # test out another general k vector [0., 0., 1.], in conventional settings
    # for Fm-3m. Bbelow are the manually generated satellite peak positions.
    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
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

    # here, we just need to create a dummy structure with the right space
    # group Fm-3m and use tools like `data2config` to expand the structure to
    # P1 symmetry and grab the atoms coordinates (stored in `atom_pos_p1`
    # below).
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

    # the hkl list file was generated by loading the dummy structure into
    # VESTA and using its powder diffraction simulator.
    hkl_refls = list()
    for i in range(6):
        for j in range(6):
            for k in range(6):
                hkl_refls.append([i, j, k])

    # as a double-check, we can load the dummy structure into the `seekpath`
    # web interface and it will tell us the right Bravais lattice type.
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
        kscope=[0., 0.5, 0., 1., 0., 1.5],
        processes=16
    )

    k_opt = k_search.kOptFinder()
    k_opt_final = k_search.kVecPrimToConv(k_opt[0])
    print("\nOptimal candidate of k vector\n===")
    for i, k_vec in enumerate(k_opt_final):
        msg = "k vector: [{:.5F}, {:.5F}, {:5F}]".format(
            k_vec[0], k_vec[1], k_vec[2]
        )
        msg += ", Indicator distance: {:.5E}".format(
            k_opt[1][i]
        )
        print(msg)

    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    # High symmetry direction.
    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    # test out another general k vector [.14, .0, .0], in conventional settings
    # for P/6mmm. Bbelow are the manually generated satellite peak positions.
    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
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

    # here, we just need to create a dummy structure with the right space
    # group P/6mmm and use tools like `data2config` to expand the structure to
    # P1 symmetry and grab the atoms coordinates (stored in `atom_pos_p1`
    # below).
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

    # the hkl list file was generated by loading the dummy structure into
    # VESTA and using its powder diffraction simulator.
    hkl_refls = list()
    for i in range(6):
        for j in range(6):
            for k in range(6):
                hkl_refls.append([i, j, k])

    # as a double-check, we can load the dummy structure into the `seekpath`
    # web interface and it will tell us the right Bravais lattice type.
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
        kscope=[0., 0.5, 0., 1., 0., 1.5],
        processes=16
    )

    k_opt = k_search.kOptFinder()
    k_opt_final = k_search.kVecPrimToConv(k_opt[0])
    print("\nOptimal candidate of k vector\n===")
    for i, k_vec in enumerate(k_opt_final):
        msg = "k vector: [{:.5F}, {:.5F}, {:5F}]".format(
            k_vec[0], k_vec[1], k_vec[2]
        )
        msg += ", Indicator distance: {:.5E}".format(
            k_opt[1][i]
        )
        print(msg)

    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    # General position
    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    # test out a general k vector of [0.33, 0.21, 1.5] in conventional setting
    # for R-3m. Below are the manually generated satellite peak positions.
    # +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    if gen_pos:
        lambda_val = 2.4109
        two_thetas = [
            14.721418097657041,
            29.173461414626605,
            29.52121838954299,
            33.94660472806508,
            39.850980755147816,
            40.459219511964385,
            42.89866262409193,
            49.92495034620636,
            51.56221863607914,
            54.73338120890784,
            56.33546485716554,
            59.68054849312238,
            60.8394556388351
        ]

        spos_gen = list()
        for two_theta in two_thetas:
            d_tmp = lambda_val / (2. * np.sin(np.deg2rad(two_theta) / 2.))
            spos_gen.append(d_tmp)

        # here, we just need to create a dummy structure with the right space
        # group R-3m and use tools like `data2config` to expand the structure
        # to P1 symmetry and grab the atoms coordinates (stored in
        # `atom_pos_p1` below).
        latt_params = [
            7.424860,
            7.424860,
            19.497026,
            90.000000,
            90.000000,
            120.000000
        ]
        latt_vec = kvs.lat_params_to_vec(latt_params)

        atoms_labels = [1, 1, 1]

        atom_pos_p1 = [
            [0.000000, 0.000000, 0.000000],
            [0.666667, 0.333333, 0.333333],
            [0.333333, 0.666667, 0.666667]
        ]

        # the hkl list file was generated by loading the dummy structure into
        # VESTA and using its powder diffraction simulator.
        hkl_refls = list()
        for i in range(6):
            for j in range(6):
                for k in range(6):
                    hkl_refls.append([i, j, k])

        # as a double-check, we can load the dummy structure into the
        # `seekpath` web interface and it will tell us the right Bravais
        # lattice type.
        brav_type = "hR"
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
            kscope=[0., 0.5, 0., 1., 0., 1.5],
            processes=16
        )

        k_opt = k_search.kOptFinder()
        k_opt_final = k_search.kVecPrimToConv(k_opt[0])
        print("\nOptimal candidate of k vector\n===")
        for i, k_vec in enumerate(k_opt_final):
            msg = "k vector: [{:.5F}, {:.5F}, {:5F}]".format(
                k_vec[0], k_vec[1], k_vec[2]
            )
            msg += ", Indicator distance: {:.5E}".format(
                k_opt[1][i]
            )
            print(msg)
