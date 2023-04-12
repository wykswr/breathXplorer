from typing import Tuple

import numpy as np
from scipy import interpolate
from scipy import signal
from sklearn.mixture import GaussianMixture


def time_with_peak(times: np.ndarray, values: np.ndarray, method: str, n_peak=1) -> Tuple[np.ndarray, float]:
    """
    Given an array of time points and an array of corresponding values, this function finds the indices and total time
    at which peaks occur, based on the specified method ('Topological' or 'Gaussian') and the number of peaks to
    consider (n_peak).

    :param times: An array of time points.
    :param values: An array of corresponding values (e.g. concentration or intensity).
    :param method: The method used to find peaks ('Topological' or 'Gaussian').
    :param n_peak: The number of peaks to consider.
    :return: A tuple containing the unique indices at peaks and the total time spent at peaks.
    """
    if method == 'Topological':
        range_l = find_peak(times, values)
    else:
        range_l = find_gaussian(times, values, n_peak)
    at_peak = list()
    total_time = 0
    for s, e, _ in range_l:
        total_time += (e - s)
        ids, ide = np.abs(times - s).argmin(), np.abs(times - e).argmin()
        for i in range(ids, ide + 1):
            at_peak.append(i)
    at_peak = np.unique(at_peak)
    return at_peak, total_time


def __lin_interp(x: np.ndarray, y: np.ndarray, i: int, cut_y: float) -> float:
    """
    Linearly interpolates the x value corresponding to a given y value (cut_y) between two points in the x-y arrays.

    :param x: An array of x values.
    :param y: An array of corresponding y values.
    :param i: The index of the first point to consider for interpolation.
    :param cut_y: The y value for which to find the corresponding x value.
    :return: The interpolated x value.
    """
    return x[i] + (x[i + 1] - x[i]) * ((cut_y - y[i]) / (y[i + 1] - y[i]))


def __find_one_range(x: np.ndarray, y: np.ndarray, idx: int) -> tuple:
    """
    Finds the range around a peak, determined by a specified index (idx).

    :param x: An array of x values.
    :param y: An array of corresponding y values.
    :param idx: The index of the peak.
    :return: A tuple containing the start and end x values of the range, and the y value at the peak.
    """
    base = min(y)
    cut_y = (y[idx] - base) * .25 + base  # can influence effect of finding range
    signs = np.sign(np.add(y, -cut_y))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    before = 0
    after = 0
    for i in zero_crossings_i:
        if i > idx:
            after = i
            break
        before = i
    return __lin_interp(x, y, before, cut_y), __lin_interp(x, y, after, cut_y), y[idx]


def __find_all_range(x: np.ndarray, y: np.ndarray) -> list:
    """
    Finds the ranges around all peaks in the x-y data.

    :param x: An array of x values.
    :param y: An array of corresponding y values.
    :return: A list of tuples containing the start and end x values of the ranges, and the y values at the peaks.
    """
    peaks = signal.argrelmax(y)[0]
    dy = y - min(y)
    peaks = peaks[dy[peaks] > max(dy) * .2]  # can influence effect of finding range
    return [__find_one_range(x, y, lb) for lb in peaks]


def __is_in(a: tuple, b: tuple) -> bool:
    """
    Determines whether two ranges overlap.

    :param a: A tuple representing the first range (start, end, value).
    :param b: A tuple representing the second range (start, end, value).
    :return: True if the ranges overlap, False otherwise.
    """
    return (a[0] <= b[0] and a[1] >= b[1]) or (a[0] >= b[0] and a[1] <= b[1])


def __merge_ranges(ori_range: list) -> list:
    """
    Merges overlapping ranges in the input list.

    :param ori_range: A list of tuples representing ranges (start, end, value).
    :return: A list of merged ranges (non-overlapping).
    """
    merged_ranges = list()
    for r in ori_range:
        merged = False
        for i, existing_range in enumerate(merged_ranges):
            if __is_in(r, existing_range):
                merged_ranges[i] = (
                    min(r[0], existing_range[0]), max(r[1], existing_range[1]), max(r[2], existing_range[2]))
                merged = True
                break
        if not merged:
            merged_ranges.append(r)
    if len(merged_ranges) == len(ori_range):
        return merged_ranges
    else:
        return __merge_ranges(merged_ranges)


def find_peak(x: np.ndarray, y: np.ndarray) -> list:
    """
    Finds the x's range of each peak using the 'Topological' method.

    :param x: An array of time points.
    :param y: An array of corresponding values (e.g. concentration or intensity).
    :return: A list of tuples containing the start and end x values of the ranges, and the y values at the peaks.
    """
    f = interpolate.interp1d(x, y, kind='cubic')
    n_x = np.union1d(np.linspace(np.min(x), np.max(x), max(1000, len(x))), x)  # can influence effect of finding range
    n_y = f(n_x)
    ori_range = __find_all_range(n_x, n_y)
    merged_ranges = __merge_ranges(ori_range)
    return merged_ranges


def find_gaussian(x: np.ndarray, y: np.ndarray, n: int) -> list:
    """
    Uses the Gaussian Mixture Model to find the range of each peak.

    :param x: An array of time points.
    :param y: An array of corresponding values (e.g. concentration or intensity).
    :param n: The number of peaks to consider.
    :return: A list of tuples containing the start and end x values of the ranges, and the y values at the peaks.
    """
    ny = y.copy()
    ny -= ny.min()
    ny /= ny.sum()
    samples = x[np.random.choice(len(ny), size=5000, p=ny)]
    gmm = GaussianMixture(n_components=n)
    gmm.fit(samples.reshape(-1, 1))
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten())
    intensities = [y[np.abs(x - m).argmin()] for m in means]
    return [(m - std, m + std, intensity) for m, std, intensity in zip(means, stds, intensities)]
