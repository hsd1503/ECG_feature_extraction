"""
特征提取器
"""
import os

import numpy as np
from tools.ResampleTools import adaptive_resampling
from dtw import accelerated_dtw
import scipy.signal as signal


def morphological_variability(clean_time_sequence=None, fs=None, peaks=None, head_weight=0.4, tail_weight=0.6,
                              segments: np.ndarray = None, target_seg_len=50, frequency_bounds: list = [0.3, 0.55],
                              base_template=None, dist='euclidean'):
    """

    Parameters
    ------------------------
    clean_time_sequence: 干净的时间序列
    fs: 时间序列采样率
    peaks: 时间序列对应的分析峰值点
    head_weight: 峰值点前的
    tail_weight:
    target_seg_len:
    frequency_bounds: 最终需要的频带范围


    Returns
    -----------------------
    时间序列在以peaks为中心附近获取到的形态变异性
    `None`表示无法获取相关形态变异性

    """
    if segments is not None:
        fragments = segments.tolist()
    else:
        # fragment
        fragments = extract_beats(clean_time_sequence, peaks, head_weight, tail_weight)
        if fragments is None:
            return None

    # record the resample time
    time_records = []
    # resampling
    for idx in range(len(fragments)):
        time_, resampling = adaptive_resampling(fragments[idx], target_seg_len)
        time_ *= (1.0 / fs)
        fragments[idx] = resampling
        time_records.append(time_)

    if base_template is not None:
        time_, resampling = adaptive_resampling(base_template, target_seg_len)
        time_ *= (1.0 / fs)
        fragments.insert(0, resampling)
        time_records.insert(0, time_)

    distances = []
    time_infos = []

    if base_template is not None:
        for idx in range(1, len(fragments)):
            d, cost_matrix, acc_cost_matrix, path = accelerated_dtw(fragments[0], fragments[idx],
                                                                    dist=dist,
                                                                    warp=3)

            costs = np.array([cost_matrix[path[0][_]][path[1][_]] for _ in range(len(path[0]))])
            distances.append(costs)
            max_ = max(max(time_records[0]), max(time_records[idx]))
            min_ = min(min(time_records[0]), min(time_records[idx]))
            time_infos.append(np.linspace(min_, max_, len(path[0])))
    else:
        for idx in range(len(fragments) - 1):
            d, cost_matrix, acc_cost_matrix, path = accelerated_dtw(fragments[idx], fragments[idx + 1],
                                                                    dist=dist,
                                                                    warp=3)
            costs = np.array([cost_matrix[path[0][_]][path[1][_]] for _ in range(len(path[0]))])
            distances.append(costs)
            max_ = max(max(time_records[idx]), max(time_records[idx + 1]))
            min_ = min(min(time_records[idx]), min(time_records[idx + 1]))
            time_infos.append(np.linspace(min_, max_, len(path[0])))

    mv = []
    freqs = np.linspace(0.01, 1, 1000)
    need_ = np.where(np.array(freqs >= frequency_bounds[0]) & np.array(freqs <= frequency_bounds[1]))
    for idx in range(len(distances)):
        # 频率谱
        pgram = signal.lombscargle(time_infos[idx], distances[idx], freqs, normalize=False)
        mv.append(np.sum(pgram[need_] ** 2 * (freqs[1] - freqs[0])))
    return mv


def extract_beats(clean_time_sequence, peaks, head_weight=0.4, tail_weight=0.6, need_beats_peak=False):
    """

    Parameters
    -------------------
    clean_time_sequence:
    fs:
    peaks:
    head_weight:
    tail_weight:



    Returns
    -------------------

    """

    # 首尾两个peak将不被截断
    if len(peaks) < 4:
        return None
    peaks = np.array(peaks)
    beat_peaks = []
    peaks.sort()
    result = []
    for p in range(1, len(peaks) - 1):
        last_ = peaks[p - 1]
        current_ = peaks[p]
        next_ = peaks[p + 1]
        head = round(head_weight * (current_ - last_))
        head_ = int(current_ - head)
        tail_ = int(current_ + round(tail_weight * (next_ - current_)))
        result.append(clean_time_sequence[head_:tail_])
        beat_peaks.append(head)
    if need_beats_peak:
        return result, beat_peaks
    else:
        return result