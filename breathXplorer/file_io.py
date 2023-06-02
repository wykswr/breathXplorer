from functools import reduce
from pathlib import Path
from typing import Tuple, List, Union

import numpy as np
import pandas as pd
from scipy import signal, sparse

from .cluster import cluster_mz, __drop_zero_near
from .container import FeatureSet, TandemMS, Source
from .find_peak import time_with_peak


def __purify(mz: np.ndarray, ints: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Filter the m/z values and intensities by removing noise.
    :param mz: m/z values
    :param ints: Intensities
    :return: Filtered m/z values and intensities
    """
    peaks = signal.argrelmax(ints)[0]
    peaks = peaks[ints[peaks] > ints.max() / 1000]
    return mz[peaks], ints[peaks]


def __scan_mz(file: str, line: bool) -> np.ndarray:
    """
    Read all m/z values at any time from an MS file (mzML or mzXML).
    :param file: MS file (mzML or mzXML) path
    :param line: Whether the spectra are line spectra
    :return: Unique m/z values
    """
    concat_func = lambda x, y: np.concatenate([x, y])
    source = Source(file)
    all_mz = reduce(concat_func, [sp.mz for sp in source if sp.level == 1]) if line else \
        reduce(concat_func, [__purify(sp.mz, sp.intensity)[0] for sp in source if sp.level == 1])

    return np.unique(all_mz)


def __scan_time(file: str) -> np.ndarray:
    """
    Read all scan times in an MS file (mzML or mzXML), usually in minutes.
    :param file: MS file (mzML or mzXML) path
    :return: Scan times
    """
    source = Source(file)
    all_time = [sp.scan_start_time for sp in source if sp.level == 1]
    return np.array(all_time)


def __scan_tic(file: str) -> np.ndarray:
    """
    Get the total ion current (TIC) at each rotation time.
    :param file: MS file (mzML or mzXML) path
    :return: TIC values
    """
    source = Source(file)
    tic = [sp.tic for sp in source if sp.level == 1]
    return np.array(tic)


def read_sparse(file: str, line: bool) -> pd.DataFrame:
    """
    Read all information of an MS file (mzML or mzXML), with rows as unique m/z values and columns as time.
    :param file: MS file (mzML or mzXML) path
    :param line: Whether the spectra are line spectra
    :return: Sparse DataFrame
    """
    times = __scan_time(file)
    times = pd.Series(index=times, data=np.arange(times.size))
    mzs = __scan_mz(file, line)
    mzs = pd.Series(index=mzs, data=np.arange(mzs.size))
    data = sparse.lil_matrix((times.size, mzs.size))

    source = Source(file)
    for sp in source:
        if sp.level == 2:
            continue
        mz_arr, ints_arr = (sp.mz, sp.intensity) if line else __purify(sp.mz, sp.intensity)
        data[times[sp.scan_start_time], mzs[mz_arr]] = ints_arr

    tb = pd.DataFrame.sparse.from_spmatrix(data.transpose())
    tb.index = mzs.index
    tb.columns = times.index
    return tb


def cluster_ms(ms: str, line: bool, dense: float, method: str, n_peak=1) -> dict:
    """
    Find all bars in the time to m/z heatmap.:param ms: MS file (mzML or mzXML) path
    :param line: Whether the MS file (mzML or mzXML) is a line spectra mzML
    :param dense: Control the quality of peak
    :param method: Peak picking method
    :param n_peak: Number of peaks
    :return: Dictionary with time, intensities, TIC, and peak time information
    """
    tb = read_sparse(ms, line)
    times = tb.columns.values
    mzs = tb.index.values
    tic = __scan_tic(ms)
    at_peak, total_time = time_with_peak(times, tic, method, n_peak)
    no_zero = np.count_nonzero(tb.iloc[:, at_peak], axis=1)
    projected_mz = reduce((lambda x, y: np.concatenate([x, y])),
                          [np.ones(no_zero[r]) * mzs[r] for r in range(mzs.size)])
    labels = cluster_mz(projected_mz, 0.001, int(at_peak.size * dense))
    new_mz, records = list(), list()
    for c in np.unique(labels[labels != -1]):
        this_mz = projected_mz[labels == c]
        new_mz.append(this_mz.mean())
        records.append(__drop_zero_near(tb.loc[np.unique(this_mz)].values, this_mz.mean(), np.unique(this_mz)))
    return {'time': times, "intensities": dict(zip(new_mz, records)), 'tic': tic, 'peak_time': total_time}


def gen_df(scanned_file: dict, new_inf: List[Tuple]) -> pd.DataFrame:
    """Generate the Dataframe.
    :param scanned_file: Dictionary read from scan_file
    :param new_inf: New information to append, such as class, mean m/z, correlation with CO2, etc.
                [tuple1(name1, dc1), tuple2(name2, dc2<key: m/z>), ...]
    :return: DataFrame to be saved in result.csv
    """
    mat = np.array(list(scanned_file['intensities'].values()))
    otb = pd.DataFrame(data=mat, index=list(scanned_file['intensities'].keys()), columns=scanned_file['time'])
    for name, ndc in new_inf:
        ntb = pd.DataFrame({name: list(ndc.values())}, index=list(ndc.keys()))
        otb = pd.concat([ntb, otb], axis='columns', ignore_index=False)
    return otb


def retrieve_tandem(file: Union[str, Path], feature: FeatureSet, radium: float) -> TandemMS:
    """
    Retrieve tandem MS information from an MS file (mzML or mzXML).

    :param file: MS file (mzML or mzXML) path
    :param feature: result FeatureSet
    :param radium: Radium of the feature
    :return: TandemMS object
    """
    tandem = TandemMS(feature.mz)
    tandem.build(Source(file), radium)
    return tandem
