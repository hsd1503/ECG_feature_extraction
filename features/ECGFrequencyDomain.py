import numpy as np
from scipy.signal import periodogram
from tools.ResampleTools import resample_interp


def psd(ecg, fs):
    """
    extract the power spectral density（PSD） of the ecg

    Parameters
    ----------
    ecg: 1d ecg serial
    fs: sampling rate of ecg

    Returns
    -------
    f: frequencies of the PSD
    psd: the PSD

    """
    fixed_fs = 300
    ts = resample_interp(ecg, fs, fixed_fs)
    f, psd = periodogram(ts, fixed_fs)
    return np.vstack([f, psd]), len(ts), fixed_fs
