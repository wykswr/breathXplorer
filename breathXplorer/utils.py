from functools import reduce
from typing import Sequence
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from numba import njit


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
    Make negative values in a numpy array 0.
    :param x: A numpy array.
    :return: A numpy array with negative values 0.
    """
    x[x < 0] = 0
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
            interp1d(original_time, original_int, kind='cubic', bounds_error=False, fill_value='extrapolate')(time))
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
