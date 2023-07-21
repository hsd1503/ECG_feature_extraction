import numpy as np


def normalize_sig_hist(sig, TH=0.01):
    """
    按信号的直方图百分比进行数据归一化

    Parameters
    ----------
    sig
    TH

    Returns
    -------

    """
    hist, edges = np.histogram(sig, 100, density=True)

    F = np.cumsum(hist)
    if F[-1] == 0:
        return sig
    F /= F[-1]
    v0 = edges[np.nonzero(F > TH)[0][0]]
    v1 = edges[np.nonzero(F < (1 - TH))[0][-1]]
    nrm = v1 - v0
    if nrm == 0:
        return sig
    sig /= float(nrm)
    return sig
