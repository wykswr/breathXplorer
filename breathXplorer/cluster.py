from functools import reduce

import numpy as np
from numba import njit
from sklearn.cluster import DBSCAN


@njit
def __drop_zero_near(arr: np.ndarray, target: float, mz: np.ndarray) -> np.ndarray:
    """
    Calculate the average of each column's nonzero element in a 2D array.

    :param arr: Input m x n matrix.
    :param target: Cluster center.
    :param mz: Input array.
    :return: 1-dimensional array with the computed values.
    """
    if arr.ndim == 1:
        return arr

    result = np.zeros(arr.shape[1])
    for i in range(arr.shape[1]):
        no_zero = arr[:, i].nonzero()[0]
        if no_zero.size == 0:
            result[i] = 0
        else:
            col_int = arr[:, i][no_zero]
            col_mz = mz[no_zero]
            result[i] = col_int[np.argmin(np.abs(col_mz - target))]
    return result


def cluster_merge(dc_ls: list) -> dict:
    """
    Merge dictionaries in dc_ls according to their keys.

    :param dc_ls: List of dictionaries with m/z as keys and values.
    :return: Dictionary with m/z as keys and arrays as values.
    """
    n_file = len(dc_ls)
    mzs = reduce((lambda x, y: np.concatenate([x, y])), map(lambda x: np.array(list(x.keys())), dc_ls))
    ints = reduce((lambda x, y: np.concatenate([x, y])), map(lambda x: np.array(list(x.values())), dc_ls))
    file_id = reduce((lambda x, y: np.concatenate([x, y])), [np.ones(len(dc_ls[i])) * i for i in range(n_file)])
    file_id = file_id.astype(int)
    clustering = DBSCAN(eps=0.0005, min_samples=2)
    labels = clustering.fit_predict(mzs.reshape(-1, 1))
    result = dict()

    for i in np.unique(labels):
        if i == -1:
            continue
        this_mzs, this_ints, this_file = map(lambda x: x[labels == i], (mzs, ints, file_id))
        mean_mz = np.mean(this_mzs)
        record = np.ones(n_file) * float('nan')
        for j in range(n_file):
            match_int = this_ints[this_file == j]
            if match_int.size == 0:
                continue
            match_mz = this_mzs[this_file == j]
            record[j] = match_int[np.argmin(np.abs(match_mz - mean_mz))]
        result[mean_mz] = record

    mzs, ints, file_id = map(lambda x: x[labels == -1], (mzs, ints, file_id))
    for m, i, f in zip(mzs, ints, file_id):
        record = np.ones(n_file) * float('nan')
        record[f] = i
        result[m] = record
    return dict(sorted(result.items(), key=lambda x: x[0]))


def cluster_mz(mzs: np.ndarray, eps: float, dense: int) -> np.ndarray:
    """
    Perform DBSCAN clustering on the input m/z array.

    :param mzs: Input m/z array.
    :param eps: The maximum distance between two samples for them to be considered as in the same neighborhood.
    :param dense: The number of samples in a neighborhood for a point to be considered as a core point.
    :return: Array of cluster labels.
    """
    clustering = DBSCAN(eps=eps, min_samples=dense)
    labels = clustering.fit_predict(mzs.reshape(-1, 1))
    return labels
