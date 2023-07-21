import numpy as np
import scipy.signal as sps


def extract_avg_wave(ecg, r_peaks, fs):
    """
    提取平均波形

    Parameters
    ----------
    ecg : ecg信号
    r_peaks : R波位置
    fs : 采样率

    Returns
    -------
    平均波形及平均波形的R波位置
    """
    b, a = sps.butter(N=4, Wn=[0.75, 35], btype='bandpass', fs=fs)
    filtered = sps.filtfilt(b, a, ecg)
    avg_len = int(np.median(np.diff(r_peaks)))
    weight = 0.35
    ahead_len = int(round(weight * avg_len))
    tail_len = int(round((1 - weight + 0.15) * avg_len))
    cache_beats = np.zeros((len(r_peaks) - 2, ahead_len+tail_len))
    for i in np.arange(1, len(r_peaks) - 1):
        start = max(r_peaks[i] - ahead_len, 0)
        end = min(r_peaks[i] + tail_len, len(ecg))
        cache_beats[i - 1][ahead_len - (r_peaks[i] - start):ahead_len] = filtered[start:r_peaks[i]]
        cache_beats[i - 1][ahead_len:ahead_len + (end - r_peaks[i])] = filtered[r_peaks[i]:end]
    return np.mean(cache_beats, axis=0), ahead_len
