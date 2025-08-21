# -*- coding: utf-8 -*-
"""k_vec_solve.py
Routines for the determination of the k-vector form compatible with ISODISTORT.

Author: Yuanpeng Zhang, with the help from Prof. Branton Campbell and Dr. Brian Toby
Date: 08/20/2025 11:56:28 EST
"""
import re
import json
import numpy as np
from fractions import Fraction
from sympy import Eq, solve, Matrix, symbols, sympify

a, b, g = symbols('a b g')


def parse_expression_string(expr_string):
    """
    Parse a string like '(a+1/2, 2b/3, g/3+2/3)' into a list of sympy expressions.

    Args:
        expr_string (str): String containing comma-separated expressions in parentheses

    Returns:
        list: List of sympy expressions
    """
    clean_string = expr_string.strip('()')
    expressions = [expr.strip() for expr in clean_string.split(',')]

    sympy_expressions = []
    for expr in expressions:
        expr_with_mult = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', expr)

        sympy_expr = sympify(expr_with_mult)
        sympy_expressions.append(sympy_expr)

    return sympy_expressions


def transform_spg_ops_to_matrix(spg_ops):
    '''Transforming space group operators in the string form into a matrix.

    Args:
        spg_ops (list): A list of operators list like ["-x + y + 1/2", "-x", "z"]

    Returns:
        numpy.ndarray: A numpy array of space group operator matrices.
    '''
    base_vectors = {
        'x': [1, 0, 0],
        'y': [0, 1, 0],
        'z': [0, 0, 1],
    }

    matrices = []

    def parse_element(element):
        resulatt_type_vector = [0, 0, 0]
        element = element.replace(' ', '')
        terms = re.split('(?=[+-])', element)

        for term in terms:
            if 'x' in term:
                sign = -1 if term.startswith('-') else 1
                resulatt_type_vector = [rv + sign * bv for rv, bv in zip(resulatt_type_vector, base_vectors['x'])]
            elif 'y' in term:
                sign = -1 if term.startswith('-') else 1
                resulatt_type_vector = [rv + sign * bv for rv, bv in zip(resulatt_type_vector, base_vectors['y'])]
            elif 'z' in term:
                sign = -1 if term.startswith('-') else 1
                resulatt_type_vector = [rv + sign * bv for rv, bv in zip(resulatt_type_vector, base_vectors['z'])]

        return resulatt_type_vector

    def extract_numeric_part(expression):
        numeric_pattern = r'[-]?\d+\/\d+|\d+\.\d+|\d+'
        numeric_parts = re.findall(numeric_pattern, expression)

        return numeric_parts[0] if numeric_parts else '0'

    for operation in spg_ops:
        matrix = []

        for element in operation:
            vector = [float(item) for item in parse_element(element)]
            num_part = float(Fraction(extract_numeric_part(element)))
            vector.append(num_part)
            matrix.append(vector)

        matrix.append([0, 0, 0, 1])
        matrix_a = np.array(matrix).T
        matrices.append(matrix_a)

    return np.array(matrices)


def grab_user_conv_transform(trans_str):
    """Given the transformation string from the space group unit cell setting to
    the conventional setting, we want to figure out the transformation matrix.

    Args:
        trans_str (str): transformation string from a CIF file with the key of '_space_group.transform_Pp_abc'

    Returns:
        numpy.ndarray: A numpy array of the transformation matrix.
    """
    x_op_1 = trans_str.split(";")[0].split(",")[0]
    x_op_2 = trans_str.split(";")[1].split(",")[0]
    if float(x_op_2) >= 0.:
        x_op_2 = f"+ {x_op_2}"
    x_op = x_op_1 + x_op_2

    y_op_1 = trans_str.split(";")[0].split(",")[1]
    y_op_2 = trans_str.split(";")[1].split(",")[1]
    if float(y_op_2) >= 0.:
        y_op_2 = f"+ {y_op_2}"
    y_op = y_op_1 + y_op_2

    z_op_1 = trans_str.split(";")[0].split(",")[2]
    z_op_2 = trans_str.split(";")[1].split(",")[2]
    if float(z_op_2) >= 0.:
        z_op_2 = f"+ {z_op_2}"
    z_op = z_op_1 + z_op_2

    op_str = [[x_op, y_op, z_op]]

    user_to_conv = transform_spg_ops_to_matrix(op_str)[0]

    return user_to_conv


def grab_user_conv_line(cif_file):
    """Grab the user to conventional setting transformation from the
    '_space_group.transform_Pp_abc' line in a given CIF file.

    Args:
        cif_file (str): Full path to the input CIF file.

    Returns:
        str: The user to conventional setting transformation string, or None if not found.
    """
    trans_str = None
    try:
        with open(cif_file, 'r', encoding='utf-8') as file:
            for line in file:
                if '_space_group.transform_Pp_abc' in line:
                    parts = line.split()
                    if len(parts) > 1:
                        trans_str = parts[1].strip()
                    break
    except FileNotFoundError:
        print(f"The file {cif_file} was not found.")

    return trans_str


def latt_prim_bases_ext(latt_type):
    """Prime lattice basis for a given lattice type.

    Args:
        latt_type (str): Lattice type

    Returns:
        list: Prime lattice basis
    """
    if latt_type == "A":
        lp_matrix = np.array([[1, 0, 0], [0, 1/2, 1/2], [0, -1/2, 1/2]])
    elif latt_type == "B":
        lp_matrix = np.array([[1/2, 0, -1/2], [0, 1, 0], [1/2, 0, 1/2]])
    elif latt_type == "C":
        lp_matrix = np.array([[1/2, 1/2, 0], [-1/2, 1/2, 0], [0, 0, 1]])
    elif latt_type == "F":
        lp_matrix = np.array([[0, 1/2, 1/2], [1/2, 0, 1/2], [1/2, 1/2, 0]])
    elif latt_type == "I":
        lp_matrix = np.array([[-1/2, 1/2, 1/2], [1/2, -1/2, 1/2], [1/2, 1/2, -1/2]])
    elif latt_type == "R":
        lp_matrix = np.array([[2/3, 1/3, 1/3], [-1/3, 1/3, 1/3], [-1/3, -2/3, 1/3]])
    else:
        lp_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    lp_matrix_t = lp_matrix.T
    lp_matrix_ext = np.hstack((lp_matrix_t, np.zeros((lp_matrix_t.shape[0], 1))))
    lp_matrix_ext = np.vstack((lp_matrix_ext, np.array([0, 0, 0, 1])))

    return lp_matrix_ext


def solve_kvec_params(numkvec, kform):
    """Given input k vector and the alternative k vector form, try to
    solve the parameters associated with the provided form.

    Args:
        numkvec (list): k vector
        kform (list): Candidate k vector form

    Returns:
        list or dict: Solution for parameters involved in the provided k vector form.
            If a solution is found, the return will be a dictionary with parameter names as keys and their
            corresponding values as dict values. If no solution is found, an empty list is returned.
    """
    kform_m = Matrix(kform)
    kvars = list(kform_m.free_symbols)
    kequats = [Eq(kform_m[i] - n, 0) for i, n in enumerate(numkvec)]
    ksols = solve(kequats, kvars)

    return ksols


def find_kvec_form_param(spg_num, isocif_cif, k_vec, k_forms):
    """Find the parameters for a given k vector form.

    Args:
        spg_num (int): Space group number.
        isocif_cif (str): Path to the CIF file exported from ISOCIF.
        k_vec (list): k vector
        k_forms (list): List of symbolic alternative k vector forms.
            Entries should be in the form of sympy expressions.

    Returns:
        tuple: (index of k_form, a, b, g)
    """
    with open("./sgtables.json", "r") as f:
        sgtables = json.load(f)

    spg_sym = sgtables[str(spg_num)]["SPGSym"]
    user_aff_ops = transform_spg_ops_to_matrix(sgtables[str(spg_num)]["SPGOps"])

    user_to_ref_conv_p_l = grab_user_conv_line(isocif_cif)
    user_to_ref_conv_p = grab_user_conv_transform(user_to_ref_conv_p_l)
    user_to_ref_conv_q = np.linalg.inv(user_to_ref_conv_p)

    conv_to_pr_p = latt_prim_bases_ext(spg_sym[0])
    conv_to_pr_q = np.linalg.inv(conv_to_pr_p)

    user_to_ref_pr_q = np.dot(conv_to_pr_q, user_to_ref_conv_q)
    user_to_ref_pr_p = np.linalg.inv(user_to_ref_pr_q)

    ref_pt_ops = np.matmul(user_to_ref_pr_q, np.matmul(user_aff_ops, user_to_ref_pr_p))

    k_vec_pr = np.dot(np.array(k_vec), user_to_ref_pr_p[0:3, 0:3])
    k_forms_pr = np.matmul(k_forms, user_to_ref_pr_p[0:3, 0:3])

    def sort_keys(dict_entry):
        # Sort keys by priority: a, then b, then g
        key_priority = {a: 1, b: 2, g: 3}
        return [key_priority.get(key, float('inf')) for key in dict_entry.keys()]

    k_vec_pr_arms = np.matmul(k_vec_pr, ref_pt_ops[0:3, 0:3])
    alt_solutions = []
    for arm in k_vec_pr_arms:
        arm_sols = []
        for k_form in k_forms_pr:
            sol = solve_kvec_params(arm, k_form)
            if sol:
                arm_sols.append(sol)
            else:
                arm_sols.append(None)

        min_keys = float('inf')
        best_dict = None
        best_index = -1

        for i, entry in enumerate(arm_sols):
            if isinstance(entry, dict):
                num_keys = len(entry)

                # Compare this entry to the current best entry
                if num_keys < min_keys or (num_keys == min_keys and sort_keys(entry) < sort_keys(best_dict)):
                    min_keys = num_keys
                    best_dict = entry
                    best_index = i

        alt_solutions.append([best_index, min_keys, best_dict])

    min_keys = float('inf')
    best_dict = None
    best_index = -1

    for entry in alt_solutions:
        if entry[0] != -1:
            num_keys = entry[1]

            if num_keys < min_keys or (num_keys == min_keys and sort_keys(entry[2]) < sort_keys(best_dict)):
                min_keys = num_keys
                best_dict = entry[2]
                best_index = entry[0]

    return (best_index, best_dict)
