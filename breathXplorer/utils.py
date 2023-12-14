from functools import reduce
from typing import Sequence

import numpy as np
import pandas as pd
from numba import njit
from scipy.interpolate import interp1d


def cal_auc(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate the area under the curve.
    :param x: x values.
    :param y: y values.
    :return: Area under the curve.
    """
    f = interp1d(x, y, kind='cubic')
    nx = np.union1d(np.linspace(np.min(x), np.max(x), max(800, len(x))), x)
    ny = f(nx)
    return np.trapz(ny, nx)


def score(scanned: dict, factor=1.0) -> dict:
    """
    Calculate the area under the curve for each m/z.
    :param scanned: Scanned data.
    :param factor: Factor to divide the area under the curve.
    :return: A dictionary of m/z and area under the curve.
    """
    times = scanned['time']
    mzs = scanned['intensities'].keys()
    return dict(zip(mzs, [cal_auc(times, v) / factor for v in scanned['intensities'].values()]))


# make negative in  np.ndarray 0
@njit
def __make_zero(x: np.ndarray) -> np.ndarray:
    """
    Make negative float values in a numpy array 0.
    :param x: A numpy array.
    :return: A numpy array with negative values 0.
    """
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = 0
    return x


def interpolate_time(tb: pd.DataFrame, time: np.ndarray) -> pd.DataFrame:
    """
    Interpolate the intensities at given time points.
    :param tb: Extracted DataFrame from one sample.
    :param time: Time points to interpolate.
    :return: Interpolated DataFrame.
    """
    intensity = tb.intensity
    tb = tb.iloc[:, 1:]
    mz = tb.index
    original_time = tb.columns.values.astype(float)
    mat = [
        __make_zero(
            interp1d(original_time, original_int, kind='linear', bounds_error=False, fill_value='extrapolate')(time))
        for original_int in tb.values]
    tb = pd.DataFrame(mat, index=mz, columns=time)
    tb['intensity'] = intensity
    # make intensity the first column
    cols = tb.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    tb = tb[cols]
    return tb


def time_union(tbs: Sequence[pd.DataFrame]) -> np.ndarray:
    """
    Get the union of time points from a list of DataFrames.
    :param tbs: List of Extracted DataFrames.
    :return: Union of time points.
    """
    return reduce(lambda x, y: np.union1d(x, y), [tb.columns[1:].values.astype(float) for tb in tbs])


_adducts = {
    'M+H': 1.007276,
    'M+H-H2O': -17.00384,
    'M+H+H2O': 19.01839,
    'M+Na': 23.98922,
}


def annotate_adduct(mz_values: np.ndarray, thres: float) -> np.ndarray:
    """
    Annotate adducts for a feature table.

    :param mz_values: m/z values.
    :param thres: Threshold for m/z difference.
    :return: Tuple of adducts and their parent ion m/z's indices.
    """

    def describe(adduct: str, idx: int) -> str:
        return f'{adduct} of m/z:{round(mz_values[idx], 4)}' if idx != -1 else 'unknown adduct'

    adducts = np.array(['unknown adduct' for _ in range(mz_values.size)])
    parent = np.ones(shape=(mz_values.size,), dtype=int) * -1

    # compare the pair of m/z values to see if they are close enough based on the adducts
    for i in range(mz_values.size):
        if parent[i] != -1:
            continue
        for j in range(i + 1, len(mz_values)):
            if parent[j] != -1:
                continue
            mz_i = mz_values[i]
            mz_j = mz_values[j]
            for adduct, mass in _adducts.items():
                if abs(mz_j - mz_i - mass) < thres:
                    adducts[j] = adduct
                    parent[j] = i
                    break

    return np.array([describe(adduct, idx) for adduct, idx in zip(adducts, parent)])


_isotopes = {
    'M+1': 1.00335,
    'M+2': 2.00671,
}


def annotate_isotope(mz_values: np.ndarray, thres: float) -> np.ndarray:
    """
    Annotate isotopes for a feature table.

    :param mz_values: m/z values.
    :param thres: Threshold for m/z difference.
    :return: Tuple of isotopes ion m/z's indices.
    """

    def describe(isotope: str, idx: int) -> str:
        return f'{isotope} of m/z:{round(mz_values[idx], 4)}' if idx != -1 else 'unknown isotope'

    isotopes = np.array(['unknown isotope' for _ in range(mz_values.size)])
    parent = np.ones(shape=(mz_values.size,), dtype=int) * -1

    # compare the pair of m/z values to see if they are close enough based on the isotopes

    for i in range(mz_values.size):
        if parent[i] != -1:
            continue
        for j in range(i + 1, len(mz_values)):
            if parent[j] != -1:
                continue
            mz_i = mz_values[i]
            mz_j = mz_values[j]
            for isotope, mass in _isotopes.items():
                if abs(mz_j - mz_i - mass) < thres:
                    isotopes[j] = isotope
                    parent[j] = i
                    break

    return np.array([describe(isotope, idx) for isotope, idx in zip(isotopes, parent)])
