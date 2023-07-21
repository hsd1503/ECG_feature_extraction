"""
HRV频域特征
"""
from scipy.signal import welch
from scipy.interpolate import interp1d
from neurokit2.signal import signal_power
import numpy as np
from features.HeartRateVariability.common import FrequencyDomainIndices, hrv_get_rri
from neurokit2 import hrv


def frequencies(peaks: list, sampling_rate=1000, rri: list = None) -> FrequencyDomainIndices:
    """
    频域分析
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `FrequencyDomainIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate, interpolate=True)
    fdis = FrequencyDomainIndices()
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]
    # resample rri
    x = np.cumsum(rrs)
    x -= x[0]
    x /= 1000
    fs = fdis.fs
    resample_x = np.arange(x[0], x[-1], 1 / fs)
    fun = interp1d(x, rrs)
    resample_rrs = fun(resample_x)
    # power spectrum
    fxx, pxx = welch(x=resample_rrs, fs=fs)
    frequency_band = [fdis.ulf_band, fdis.vlf_band, fdis.lf_band, fdis.hf_band, fdis.vhf_band]
    ulf_indexes = np.logical_and(fxx >= frequency_band[0][0], fxx < frequency_band[0][1])
    fdis.ulf_power = np.trapz(y=pxx[ulf_indexes], x=fxx[ulf_indexes])
    vlf_indexes = np.logical_and(fxx >= frequency_band[1][0], fxx < frequency_band[1][1])
    fdis.vlf_power = np.trapz(y=pxx[vlf_indexes], x=fxx[vlf_indexes])
    lf_indexes = np.logical_and(fxx >= frequency_band[2][0], fxx < frequency_band[2][1])
    fdis.lf_power = np.trapz(y=pxx[lf_indexes], x=fxx[lf_indexes])
    hf_indexes = np.logical_and(fxx >= frequency_band[3][0], fxx < frequency_band[3][1])
    fdis.hf_power = np.trapz(y=pxx[hf_indexes], x=fxx[hf_indexes])
    vhf_indexes = np.logical_and(fxx >= frequency_band[4][0], fxx < frequency_band[4][1])
    fdis.vhf_power = np.trapz(y=pxx[vhf_indexes], x=fxx[vhf_indexes])
    # total power
    fdis.total_power = np.nansum([fdis.ulf_power, fdis.vlf_power, fdis.lf_power, fdis.hf_power, fdis.vhf_power])
    fdis.lf_n = fdis.lf_power / fdis.total_power
    if np.isnan(fdis.lf_n):
        fdis.lf_n = 0
    fdis.hf_n = fdis.hf_power / fdis.total_power
    if np.isnan(fdis.hf_n):
        fdis.hf_n = 0

    if fdis.hf_power == 0 and fdis.lf_power == 0:
        fdis.lf_hf = 1
    else:
        fdis.lf_hf = fdis.lf_power / fdis.hf_power

    return fdis
