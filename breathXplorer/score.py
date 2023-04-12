import numpy as np
from scipy import interpolate


def cal_auc(x: np.ndarray, y: np.ndarray) -> float:
    f = interpolate.interp1d(x, y, kind='cubic')
    nx = np.union1d(np.linspace(np.min(x), np.max(x), max(800, len(x))), x)
    ny = f(nx)
    return np.trapz(ny, nx)


def score(scanned: dict, factor=1.0) -> dict:
    times = scanned['time']
    mzs = scanned['intensities'].keys()
    return dict(zip(mzs, [cal_auc(times, v) / factor for v in scanned['intensities'].values()]))
