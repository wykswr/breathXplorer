from pathlib import Path
from typing import Sequence, List, Union

import pandas as pd

from .cluster import cluster_merge
from .file_io import cluster_ms, gen_df
from .utils import score, time_union, interpolate_time


# a single task
def find_feature(ms: Union[Path, str], line: bool, quantity: float, method: str, n_peak: int = 1) -> pd.DataFrame:
    """
    Calculate the feature table of a single mzML file.
    :param ms:  Path of the mzML file.
    :param line:  Whether to use line mode.
    :param quantity: control the quality of peak
    :param method:  The method used to find peaks ('Topological' or 'Gaussian').
    :param n_peak:  Number of peaks to be picked
    :return:  Feature table
    """
    ms = Path(ms)
    scanned = cluster_ms(str(ms.absolute()), line, quantity, method, n_peak)
    try:
        scores = score(scanned, scanned['peak_time'])
    except ZeroDivisionError:
        scores = score(scanned)
    return gen_df(scanned, [('intensity', scores)])


# merge all the result files of single tasks in the target folder
def merge_result(tbs: Sequence[pd.DataFrame], names: List[str]) -> pd.DataFrame:
    """
    Merge the feature tables of multiple mzML files.
    :param tbs:  Feature tables
    :param names:  Names of the mzML files
    :return:  Merged feature table
    """
    sub_results = [dict(zip(tb.index, tb['intensity'])) for tb in tbs]
    result = cluster_merge(sub_results)
    return pd.DataFrame(data=list(result.values()), index=list(result.keys()), columns=names)


def time_align(tbs: Sequence[pd.DataFrame]) -> Sequence[pd.DataFrame]:
    """
    Align the time of multiple feature tables.
    :param tbs:  Feature tables
    :return:  Aligned feature tables
    """
    if len(tbs) == 1:
        return tbs
    common_time = time_union(tbs)
    return [interpolate_time(tb, common_time) for tb in tbs]
