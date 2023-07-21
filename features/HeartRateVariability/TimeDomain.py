"""
时域特征
"""
import numpy as np
from scipy import stats
from features.HeartRateVariability.common import TimeDomainBasicIndices, TimeDomainStandardIndices, hrv_get_rri
from features.HeartRateVariability.time_domain import nnXX


def basic_statistics(peaks: list, sampling_rate=1000, rri: list = None) -> TimeDomainBasicIndices:
    """
    RR间隙基础统计结果
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `TimeDomainBasicIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    tdbis = TimeDomainBasicIndices()
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]
    tdbis.max_rr = np.max(rrs)
    tdbis.min_rr = np.min(rrs)
    tdbis.mean_rr = np.mean(rrs)
    tdbis.median_rr = np.median(rrs)
    tdbis.mad_rr = _mad(rrs)
    tdbis.mcv_rr = tdbis.mad_rr / tdbis.median_rr
    tdbis.count_rr = len(rrs)
    tdbis.variance_rr = np.var(rrs)
    tdbis.std_rr = np.std(rrs)
    tdbis.amplitude_rr = tdbis.max_rr - tdbis.min_rr
    percentiles = np.percentile(rrs, [5, 25, 75, 95])
    tdbis.percentile_5 = percentiles[0]
    tdbis.percentile_25 = percentiles[1]
    tdbis.percentile_75 = percentiles[2]
    tdbis.percentile_95 = percentiles[3]
    tdbis.percentile_range_75_25 = tdbis.percentile_75 - tdbis.percentile_25
    tdbis.percentile_range_95_5 = tdbis.percentile_95 - tdbis.percentile_5
    tdbis.skew_rr = stats.skew(rrs)
    tdbis.kurtosis_rr = stats.kurtosis(rrs)
    return tdbis


def hrv(peaks: list, sampling_rate=1000, rri: list = None) -> TimeDomainStandardIndices:
    """
    HRV统计结果
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `TimeDomainStandardIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]
    mean_rr = np.mean(rrs)
    sis = TimeDomainStandardIndices()
    sis.sdnn = np.std(rrs, ddof=1)
    sis.sdann_1 = _sdann(rrs, window=1)
    sis.sdnni_1 = _sdnni(rrs, window=1)
    sis.sdann_2 = _sdann(rrs, window=2)
    sis.sdnni_2 = _sdnni(rrs, window=2)
    sis.sdann_5 = _sdann(rrs, window=5)
    sis.sdnni_5 = _sdnni(rrs, window=5)
    diff_rrs = np.diff(rrs)
    sis.rmssd = np.sqrt(np.mean(diff_rrs ** 2))
    sis.sdsd = np.nanstd(diff_rrs, ddof=1)
    sis.cvnn = sis.sdnn / mean_rr
    sis.cvsd = sis.rmssd / mean_rr
    sis.pNN20 = nnXX(rrs, threshold=20)[1]
    sis.pNN50 = nnXX(rrs, threshold=50)[1]
    binsize = (1 / 128) * 1000
    bins = np.arange(0, np.max(rrs) + binsize, binsize)
    bar_y, bar_x = np.histogram(rrs, bins=bins)
    sis.tri_index = len(rrs) / np.max(bar_y)
    sis.tinn_m, sis.tinn_n = _tinn(rrs, bar_x, bar_y, binsize)
    sis.tinn = sis.tinn_m - sis.tinn_n
    return sis


def _sdann(rri, window=1):
    """
    将序列按指定时长分段，然后计算每段的均值，最后计算这些均值的标准差
    :param rri:
    :param window: 每个片段时长，单位:min
    :return:
    """
    window_size = window * 60 * 1000
    n_windows = int(np.round(np.cumsum(rri)[-1] / window_size))
    if n_windows < 3:
        return 0
    rri_cumsum = np.cumsum(rri)
    avg_rri = []
    for i in range(n_windows):
        start = i * window_size
        start_idx = np.where(rri_cumsum >= start)[0][0]
        end_idx = np.where(rri_cumsum < start + window_size)[0][-1]
        avg_rri.append(np.mean(rri[start_idx:end_idx]))
    sdann = np.nanstd(avg_rri, ddof=1)
    return sdann


def _sdnni(rri, window=1):
    """
    将序列按指定时长分段，然后计算每段的标准差，最后计算这些标准差的均值
    :param rri:
    :param window:
    :return:
    """
    window_size = window * 60 * 1000
    n_windows = int(np.round(np.cumsum(rri)[-1] / window_size))
    if n_windows < 3:
        return np.nan
    rri_cumsum = np.cumsum(rri)
    sdnn_ = []
    for i in range(n_windows):
        start = i * window_size
        start_idx = np.where(rri_cumsum >= start)[0][0]
        end_idx = np.where(rri_cumsum < start + window_size)[0][-1]
        sdnn_.append(np.nanstd(rri[start_idx:end_idx], ddof=1))
    sdnni = np.nanmean(sdnn_)
    return sdnni


def _mad(x, constant=1.4826, **kwargs):
    """Median Absolute Deviation: a "robust" version of standard deviation.

    Parameters
    ----------
    x : Union[list, np.array, pd.Series]
        A vector of values.
    constant : float
        Scale factor. Use 1.4826 for results similar to default R.

    Returns
    ----------
    float
        The MAD.

    Examples
    ----------
    >>> import neurokit2 as nk
    >>> nk.mad([2, 8, 7, 5, 4, 12, 5, 1])
    3.7064999999999997

    References
    -----------
    - https://en.wikipedia.org/wiki/Median_absolute_deviation

    """
    median = np.nanmedian(np.ma.array(x).compressed(), **kwargs)
    mad_value = np.nanmedian(np.abs(x - median), **kwargs)
    mad_value = mad_value * constant
    return mad_value


def _tinn(rri, bar_x, bar_y, binsize):
    # set pre-defined conditions
    min_error = 2 ** 14
    X = bar_x[np.argmax(bar_y)]  # bin where Y is max
    Y = np.max(bar_y)  # max value of Y
    tmp = np.where(bar_x - np.min(rri) > 0)[0]
    if len(tmp) > 0:
        n = bar_x[tmp[0]]  # starting search of N
    else:
        n = bar_x[0]
    m = X + binsize  # starting search value of M
    N = 0
    M = 0
    # start to find best values of M and N where least square is minimized
    while n < X:
        while m < np.max(rri):
            n_start = np.where(bar_x == n)[0][0]
            n_end = np.where(bar_x == X)[0][0]
            qn = np.polyval(np.polyfit([n, X], [0, Y], deg=1), bar_x[n_start:n_end + 1])
            m_start = np.where(bar_x == X)[0][0]
            m_end = np.where(bar_x == m)[0][0]
            qm = np.polyval(np.polyfit([X, m], [Y, 0], deg=1), bar_x[m_start:m_end + 1])
            q = np.zeros(len(bar_x))
            q[n_start:n_end + 1] = qn
            q[m_start:m_end + 1] = qm
            # least squares error
            error = np.sum((bar_y[n_start:m_end + 1] - q[n_start:m_end + 1]) ** 2)
            if error < min_error:
                N = n
                M = m
                min_error = error
            m += binsize
        n += binsize
    return M, N
