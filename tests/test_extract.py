from breathXplorer import cal_one


def test_cal_one(sample):
    assert abs(cal_one(sample, False, .2, "Topological", 6).shape[0] - 279) < 10
