from pathlib import Path
from typing import Sequence

import pandas as pd

from .cluster import cluster_merge
from .file_io import cluster_ms, gen_csv
from .score import score


# a single task
def cal_one(ms: Path, line: bool, quantity: float, method: str, n_peak=1) -> pd.DataFrame:
    scanned = cluster_ms(str(ms.absolute()), line, quantity, method, n_peak)
    try:
        scores = score(scanned, scanned['peak_time'])
    except ZeroDivisionError:
        scores = score(scanned)
    return gen_csv(scanned, [('intensity', scores)])


# merge all the result files of single tasks in the target folder
def merge_result(tbs: Sequence[pd.DataFrame], names: Sequence[str]) -> pd.DataFrame:
    sub_results = [dict(zip(tb.index, tb['intensity'])) for tb in tbs]
    result = cluster_merge(sub_results)
    return pd.DataFrame(data=list(result.values()), index=list(result.keys()), columns=names)
