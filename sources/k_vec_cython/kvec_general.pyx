import numpy as np
cimport numpy as np
from scipy.optimize import linear_sum_assignment
from multiprocessing import Pool
from libc.math cimport sqrt


def unique_closest(
        np.ndarray[np.float64_t, ndim=1] list1,
        np.ndarray[np.float64_t, ndim=1] list2
    ):
    cdef int i, j
    cdef np.ndarray[np.float64_t, ndim=2] cost_matrix, cost_matrix_dd
    cdef np.ndarray[np.intp_t, ndim=1] row_ind, col_ind

    cost_matrix = np.zeros((len(list2), len(list1)), dtype=np.float64)
    cost_matrix_dd = np.zeros((len(list2), len(list1)), dtype=np.float64)
    for i in range(len(list2)):
        for j in range(len(list1)):
            cost_matrix[i, j] = ((list1[j] - list2[i]) / list2[i])**2
            cost_matrix_dd[i, j] = (list1[j] - list2[i])**2

    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    mapping = {list2[i]: list1[j] for i, j in zip(row_ind, col_ind)}
    sum_of_squares = sqrt(np.sum(cost_matrix[row_ind, col_ind])) / len(list2)
    ave_dd = sqrt(np.sum(cost_matrix_dd[row_ind, col_ind])) / len(list2)
    max_dd = sqrt(np.max(cost_matrix_dd[row_ind, col_ind]))

    return (mapping, sum_of_squares, ave_dd, max_dd)


def indDistCalc(
        list point,
        list hkl_refls,
        list superpeaks,
        list rep_prim_latt
    ):
    cdef int i
    cdef np.ndarray[np.float64_t, ndim=1] satellite_peaks_p, satellite_peaks_m
    cdef np.ndarray[np.float64_t, ndim=1] satellite_peaks

    satellite_peaks_p = np.empty(len(hkl_refls), dtype=np.float64)
    satellite_peaks_m = np.empty(len(hkl_refls), dtype=np.float64)

    def calcSatelliteP(list hkl):
        cdef np.ndarray[np.float64_t, ndim=1] hkl_p_k
        cdef np.ndarray[np.float64_t, ndim=1] k_cart
        cdef np.float64_t d_hkl_p_k

        hkl_p_k = np.array(hkl) + np.array(point)
        k_cart = np.matmul(hkl_p_k, np.array(rep_prim_latt))
        d_hkl_p_k = 2. * np.pi / np.linalg.norm(k_cart)
        return d_hkl_p_k

    def calcSatelliteM(list hkl):
        cdef np.ndarray[np.float64_t, ndim=1] hkl_m_k
        cdef np.ndarray[np.float64_t, ndim=1] k_cart
        cdef np.float64_t d_hkl_m_k

        hkl_m_k = np.array(hkl) - np.array(point)
        k_cart = np.matmul(hkl_m_k, np.array(rep_prim_latt))
        d_hkl_m_k = 2. * np.pi / np.linalg.norm(k_cart)
        return d_hkl_m_k

    for i in range(len(hkl_refls)):
        satellite_peaks_p[i] = calcSatelliteP(hkl_refls[i])
        satellite_peaks_m[i] = calcSatelliteM(hkl_refls[i])

    satellite_peaks = np.concatenate(
        (satellite_peaks_p, satellite_peaks_m)
    )

    _, indicator_dist, ave_dd, max_dd = unique_closest(
        satellite_peaks,
        np.array(superpeaks, dtype=np.float64)
    )

    return (point, indicator_dist, ave_dd, max_dd)

def parallel_proc(
        np.ndarray[np.float64_t, ndim=2] points,
        list hkl_refls,
        list superpeaks,
        list rep_prim_latt,
        int processes=4
    ):
    cdef int i
    with Pool(processes=processes) as pool:
        result = pool.starmap(
            indDistCalc,
            [
                (
                    [float(point) for point in point_list],
                    hkl_refls,
                    superpeaks,
                    rep_prim_latt
                )
                for point_list in points
            ]
        )

    return sorted(result, key=lambda x: x[1])[:10]
