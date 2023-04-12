import sys
from pathlib import Path

sys.path.append(Path(__file__).parent.parent.as_posix())

from breathXplorer import find_feature


def test_cal_one(sample):
    assert abs(find_feature(sample, False, .2, "Topological", 6).shape[0] - 279) < 10
    assert abs(find_feature(sample, False, .2, "Gaussian", 6).shape[0] - 279) < 10
