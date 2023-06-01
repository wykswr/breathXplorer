import sys
from pathlib import Path

sys.path.append(Path(__file__).parent.parent.as_posix())

from breathXplorer import find_feature, merge_result


def test_cal_one(sample):
    assert abs(len(find_feature(sample, False, .2, "Topological", 6)) - 279) < 10
    assert abs(len(find_feature(sample, False, .2, "Gaussian", 6)) - 279) < 10


def test_cal_all(sample):
    tbs = [find_feature(f, False, .2, "Topological", 6) for f in [sample, sample, sample]]
    assert abs(len(merge_result(tbs, ["a", "b", "c"])) - 284) < 10
