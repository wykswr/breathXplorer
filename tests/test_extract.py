import sys
from pathlib import Path

sys.path.append(Path(__file__).parent.parent.as_posix())

from breathXplorer import find_feature, merge_result, time_align


def test_cal_one(sample):
    assert abs(find_feature(sample, False, .2, "Topological", 6).shape[0] - 279) < 10
    assert abs(find_feature(sample, False, .2, "Gaussian", 6).shape[0] - 279) < 10


def test_cal_all(sample):
    tbs = [find_feature(f, False, .2, "Topological", 6) for f in [sample, sample, sample]]
    assert abs(merge_result(tbs, ["a", "b", "c"]).shape[0] - 284) < 10


def test_time_align(sample):
    tb = find_feature(sample, False, .2, "Topological", 6)
    tb2 = tb.copy()
    tb2.columns = ["intensity"] + list(tb.columns[1:].values + 0.001)
    tb3 = tb.copy()
    tb3.columns = ["intensity"] + list(tb.columns[1:].values - 0.002)
    tbs = time_align([tb, tb2, tb3])
    assert (tbs[0].columns[1:].values - tbs[1].columns[1:].values).sum() < 0.001
