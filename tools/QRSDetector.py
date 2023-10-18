import json
import requests
import numpy as np
import neurokit2 as nk
import scipy.signal as sps
from tools.ResampleTools import resample_interp


def detect_qrs(sig: list, fs: int, url: str) -> list:
    """
    detect qrs complex index using (HTTP) API

    Parameters
    ----------
    sig : 1-D signal
    fs : signal sample rate
    url: http url

    Returns
    -------
    qrs complex index

    """
    request_body = {
        'ecgData': sig,
        'sampleRate': fs
    }
    data_json = json.dumps(request_body)
    reps = requests.post(url, data=data_json)
    # try:
    reps = json.loads(reps.text)
    # except json.decoder.JSONDecodeError as e:
    #     print()
    #     raise e
    qrs = reps.get('data', None)
    if qrs is None:
        return None
    else:
        b, a = sps.butter(N=2, Wn=[5, 30], btype='bandpass', fs=fs)
        filtered = sps.filtfilt(b, a, sig)
        fixed = []
        for _ in qrs:
            start = max(0, int(_ - 0.1 * fs))
            end = min(len(sig), int(_ + 0.1 * fs))
            fixed.append(start + np.argmax(filtered[start:end]))
        return fixed


def simple_qrs_detector(ecg, fs):
    """
    简易的QRS波群检测器

    Parameters
    ----------
    ecg : 1D ECG数据
    fs : 信号采样率

    Returns
    -------
    QRS波群位置

    """
    fixed_fs = 500
    # remove power-line interference
    b, a = sps.iirnotch(50, 50, fs)
    ecg = sps.filtfilt(b, a, ecg)
    whole_fixed_fs_ecg = resample_interp(ecg, fs, fixed_fs)
    # simple filter
    b, a = sps.butter(N=4, Wn=[0.5, 45], btype='bandpass', fs=fixed_fs)
    whole_fixed_fs_ecg = sps.filtfilt(b, a, whole_fixed_fs_ecg)
    # clean ecg
    # clean ecg
    b, a = sps.butter(N=4, Wn=[3, 40], btype='bandpass', fs=fixed_fs)
    whole_fixed_fs_filtered = sps.filtfilt(b, a, whole_fixed_fs_ecg)
    # ecg r peaks
    _, info = nk.ecg_peaks(whole_fixed_fs_filtered, fixed_fs,
                        #    smoothwindow=0.2,
                        #    avgwindow=0.9,
                        #    gradthreshweight=1.2,
                        #    minlenweight=0.3,
                        #    mindelay=0.2
                           )
    whole_fixed_fs_r_peaks = info['ECG_R_Peaks']
    if len(whole_fixed_fs_r_peaks) < 1:
        return whole_fixed_fs_r_peaks
    else:
        return [int(_ / fixed_fs * fs) for _ in whole_fixed_fs_r_peaks]
