"""
单个个Beat波形起止定位
"""
import math
import scipy
import numpy as np
from tools.ResampleTools import resample_interp
from neurokit2.signal import signal_findpeaks, signal_smooth


def signal_zerocrossings(signal, direction="both"):
    """Locate the indices where the signal crosses zero.

    Note that when the signal crosses zero between two points, the first index is returned.

    Parameters
    ----------
    signal : Union[list, np.array, pd.Series]
        The signal (i.e., a time series) in the form of a vector of values.
    direction : str
        Direction in which the signal crosses zero, can be "positive", "negative" or "both" (default).

    Returns
    -------
    array
        Vector containing the indices of zero crossings.

    Examples
    --------

    """
    df = np.diff(np.sign(signal))
    if direction in ["positive", "up"]:
        zerocrossings = np.where(df > 0)[0]
    elif direction in ["negative", "down"]:
        zerocrossings = np.where(df < 0)[0]
    else:
        zerocrossings = np.nonzero(np.abs(df) > 0)[0]

    return zerocrossings


def _dwt_adjust_parameters(sampling_rate, average_rate, duration=None, target=None):
    if target == "degree":
        # adjust defree of dwt by sampling_rate and HR
        scale_factor = (sampling_rate / 250) / (average_rate / 60)
        return int(np.log2(scale_factor))
    elif target == "duration":
        # adjust duration of search by HR
        return np.round(duration * (60 / average_rate), 3)


def _dwt_delineate_tp_peaks(heartbeat, rpeak, dwtmatr, average_rate, sampling_rate=250, qrs_width=0.13,
                            p2r_duration=0.3, rt_duration=0.3, degree_tpeak=3, degree_ppeak=2, epsilon_T_weight=0.25,
                            epsilon_P_weight=0.02):
    """
    Parameters
    ----------
    heartbeat : Union[list, np.array, pd.Series]
        The cleaned ECG channel as returned by `ecg_clean()`.
    rpeak : Union[list, np.array, pd.Series]
        The samples at which R-peaks occur. Accessible with the key "ECG_R_Peaks" in the info dictionary
        returned by `ecg_findpeaks()`.
    dwtmatr : np.array
        Output of `_dwt_compute_multiscales()`. Multiscales of wavelet transform.
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
    qrs_width : int
        Approximate duration of qrs in seconds. Default to 0.13 seconds.
    p2r_duration : int
        Approximate duration from P peaks to R peaks in seconds.
    rt_duration : int
        Approximate duration from R peaks to T peaks in secons.
    degree_tpeak : int
        Wavelet transform of scales 2**3.
    degree_tpeak : int
        Wavelet transform of scales 2**2.
    epsilon_T_weight : int
        Epsilon of RMS value of wavelet transform. Appendix (A.3).
    epsilon_P_weight : int
        Epsilon of RMS value of wavelet transform. Appendix (A.4).
    """

    srch_bndry = int(0.5 * qrs_width * sampling_rate)
    degree_add = _dwt_adjust_parameters(sampling_rate, average_rate, target="degree")
    # sanitize search duration by HR
    p2r_duration = _dwt_adjust_parameters(sampling_rate, average_rate, duration=p2r_duration, target="duration")
    rt_duration = _dwt_adjust_parameters(sampling_rate, average_rate, duration=rt_duration, target="duration")

    tpeak = math.nan
    if not np.isnan(rpeak):
        # search for T peaks from R peaks
        srch_idx_start = rpeak + srch_bndry
        srch_idx_end = min(rpeak + 2 * int(rt_duration * sampling_rate), dwtmatr.shape[1])
        dwt_local = dwtmatr[degree_tpeak + degree_add, srch_idx_start:srch_idx_end]
        height = epsilon_T_weight * np.sqrt(np.mean(np.square(dwt_local)))

        if len(dwt_local) != 0:
            ecg_local = heartbeat[srch_idx_start:srch_idx_end]
            peaks, __ = scipy.signal.find_peaks(np.abs(dwt_local), height=height)
            peaks = list(
                filter(lambda p: np.abs(dwt_local[p]) > 0.025 * max(dwt_local), peaks)
            )  # pylint: disable=W0640
            if dwt_local[0] > 0:  # just append
                peaks = [0] + peaks

            # detect morphology
            candidate_peaks = []
            candidate_peaks_scores = []
            for idx_peak, idx_peak_nxt in zip(peaks[:-1], peaks[1:]):
                correct_sign = (
                        dwt_local[idx_peak] * dwt_local[idx_peak_nxt] < 0
                )  # pylint: disable=R1716
                if correct_sign:
                    idx_zero = (
                            signal_zerocrossings(dwt_local[idx_peak: idx_peak_nxt + 1])[0] + idx_peak
                    )
                    # This is the score assigned to each peak. The peak with the highest score will be
                    # selected.
                    # score = ecg_local[idx_zero] - (float(idx_zero) / sampling_rate - (rt_duration - 0.5 * qrs_width))
                    score = abs(ecg_local[idx_zero] - np.average(heartbeat))
                    candidate_peaks.append(idx_zero)
                    candidate_peaks_scores.append(score)

            if not candidate_peaks:
                tpeak = math.nan
            else:
                tpeak = candidate_peaks[np.argmax(candidate_peaks_scores)] + srch_idx_start

    ppeak = math.nan
    if not np.isnan(rpeak):
        # search for P peaks from Rpeaks
        srch_idx_start = max(rpeak - 2 * int(p2r_duration * sampling_rate), 0)
        if srch_idx_start <= 0:
            srch_idx_start = max(rpeak - 1.5 * int(p2r_duration * sampling_rate), 0)
        if srch_idx_start <= 0:
            srch_idx_start = max(rpeak - int(p2r_duration * sampling_rate), 0)
        srch_idx_end = rpeak - srch_bndry
        dwt_local = dwtmatr[degree_ppeak + degree_add, srch_idx_start:srch_idx_end]
        height = epsilon_P_weight * np.sqrt(np.mean(np.square(dwt_local)))

        if len(dwt_local) != 0:
            ecg_local = heartbeat[srch_idx_start:srch_idx_end]
            peaks, __ = scipy.signal.find_peaks(np.abs(dwt_local), height=height)
            peaks = list(filter(lambda p: np.abs(dwt_local[p]) > 0.025 * max(dwt_local), peaks))
            if dwt_local[0] > 0:  # just append
                peaks = [0] + peaks

            # detect morphology
            candidate_peaks = []
            candidate_peaks_scores = []
            for idx_peak, idx_peak_nxt in zip(peaks[:-1], peaks[1:]):
                correct_sign = (
                        dwt_local[idx_peak] * dwt_local[idx_peak_nxt] < 0
                )  # pylint: disable=R1716
                if correct_sign:
                    idx_zero = (
                            signal_zerocrossings(dwt_local[idx_peak: idx_peak_nxt + 1])[0] + idx_peak
                    )
                    # This is the score assigned to each peak. The peak with the highest score will be
                    # selected.
                    # score = ecg_local[idx_zero] - abs(float(idx_zero) / sampling_rate - p2r_duration)
                    score = abs(ecg_local[idx_zero] - np.average(heartbeat))
                    candidate_peaks.append(idx_zero)
                    candidate_peaks_scores.append(score)

            if not candidate_peaks:
                ppeak = math.nan
            else:
                ppeak = candidate_peaks[np.argmax(candidate_peaks_scores)] + srch_idx_start

    return tpeak, ppeak


def _dwt_delineate_tp_onsets_offsets(peak, dwtmatr, average_rate, bounds, sampling_rate=250, duration_onset=0.3,
                                     duration_offset=0.3, onset_weight=0.4, offset_weight=0.4, degree_onset=2,
                                     degree_offset=2):
    onset = math.nan
    offset = math.nan
    if peak is math.nan:
        return onset, offset
    if bounds[0] <= peak <= bounds[1]:
        return onset, offset
    elif peak < bounds[0]:
        bound = bounds[0]
    else:
        bound = bounds[1]
    # sanitize search duration by HR
    duration_onset = _dwt_adjust_parameters(sampling_rate, average_rate, duration=duration_onset, target="duration")
    duration_offset = _dwt_adjust_parameters(sampling_rate, average_rate, duration=duration_offset, target="duration")
    degree = _dwt_adjust_parameters(sampling_rate, average_rate, target="degree")

    # look for onsets
    if peak < bounds[0]:
        srch_idx_start = max(peak - int(duration_onset * sampling_rate), 0)
    else:
        srch_idx_start = max(peak - int(duration_onset * sampling_rate), bound)
    srch_idx_end = peak
    if srch_idx_start >= srch_idx_end:
        pass
    else:
        dwt_local = dwtmatr[degree_onset + degree, srch_idx_start:srch_idx_end]
        onset_slope_peaks, __ = scipy.signal.find_peaks(dwt_local)
        if len(onset_slope_peaks) == 0:
            pass
        else:
            epsilon_onset = onset_weight * dwt_local[onset_slope_peaks[-1]]
            if not (dwt_local[: onset_slope_peaks[-1]] < epsilon_onset).any():
                pass
            else:
                candidate_onsets = np.where(dwt_local[: onset_slope_peaks[-1]] < epsilon_onset)[0]
                onset = candidate_onsets[-1] + srch_idx_start

    # look for offset
    srch_idx_start = peak
    if peak < bounds[0]:
        srch_idx_end = max(peak + int(duration_offset * sampling_rate), bound)
    else:
        srch_idx_end = min(max(peak + int(duration_offset * sampling_rate), bound), dwtmatr.shape[1])

    if srch_idx_start >= srch_idx_end:
        pass
    else:
        dwt_local = dwtmatr[degree_offset + degree, srch_idx_start:srch_idx_end]
        offset_slope_peaks, __ = scipy.signal.find_peaks(-dwt_local)
        if len(offset_slope_peaks) == 0:
            pass
        else:
            epsilon_offset = -offset_weight * dwt_local[offset_slope_peaks[0]]
            if not (-dwt_local[offset_slope_peaks[0]:] < epsilon_offset).any():
                pass
            else:
                candidate_offsets = (
                        np.where(-dwt_local[offset_slope_peaks[0]:] < epsilon_offset)[0]
                        + offset_slope_peaks[0]
                )
                offset = candidate_offsets[0] + srch_idx_start

    return onset, offset


def _dwt_delineate_qrs_bounds(rpeak, dwtmatr, ppeak, tpeak, average_rate, Q, S, sampling_rate=250):
    """
    定位QRS波群的相关位置

    Parameters
    ----------
    rpeak
    dwtmatr
    ppeak
    tpeak
    average_rate
    sampling_rate

    Returns
    -------

    """
    degree = _dwt_adjust_parameters(sampling_rate, average_rate, target="degree")
    onset = math.nan
    # look for onsets
    srch_idx_start = max(int(rpeak - 0.1 * sampling_rate) if ppeak is math.nan else ppeak, 0)
    srch_idx_end = rpeak
    if srch_idx_start is math.nan or srch_idx_end is math.nan:
        pass
    else:
        dwt_local = dwtmatr[2 + degree, srch_idx_start:srch_idx_end]
        onset_slope_peaks, __ = scipy.signal.find_peaks(-dwt_local)
        if len(onset_slope_peaks) == 0:
            onset = Q
        else:
            epsilon_onset = 0.5 * -dwt_local[onset_slope_peaks[-1]]
            if not (-dwt_local[: onset_slope_peaks[-1]] < epsilon_onset).any():
                onset = Q
            else:
                candidate_onsets = np.where(-dwt_local[: onset_slope_peaks[-1]] < epsilon_onset)[0]
                onset = candidate_onsets[-1] + srch_idx_start
            tmp = dwtmatr[2 + degree, int(Q - 0.02 * sampling_rate):int(Q + 0.02 * sampling_rate)]
            if len(tmp) > 3:
                if np.max(tmp) - np.min(tmp) < 0.5:
                    onset = Q

    offset = math.nan
    # look for offsets
    srch_idx_start = rpeak
    srch_idx_end = int(rpeak + 0.1 * sampling_rate) if tpeak is math.nan else tpeak
    if srch_idx_start is math.nan or srch_idx_end is math.nan:
        pass
    else:
        dwt_local = dwtmatr[2 + degree, srch_idx_start:srch_idx_end]
        onset_slope_peaks, __ = scipy.signal.find_peaks(dwt_local)
        if len(onset_slope_peaks) == 0:
            offset = S
        else:
            epsilon_offset = 0.5 * dwt_local[onset_slope_peaks[0]]
            if not (dwt_local[onset_slope_peaks[0]:] < epsilon_offset).any():
                offset = S
            else:
                candidate_offsets = (
                        np.where(dwt_local[onset_slope_peaks[0]:] < epsilon_offset)[0] + onset_slope_peaks[0]
                )
                offset = candidate_offsets[0] + srch_idx_start

    return onset, offset


def _ecg_delineator_peak_Q(rpeak, heartbeat, R, fs):
    """
    定位Q波位置

    Parameters
    ----------
    rpeak
    heartbeat
    R

    Returns
    -------

    """
    segment = np.diff(heartbeat[:rpeak])
    Q = signal_findpeaks(-1 * segment)
    if len(Q["Peaks"]) == 0:
        return math.nan, None
    Q = np.max(Q["Peaks"])
    # tmp = heartbeat[int(Q - 0.015 * fs):int(Q + 0.015 * fs)]
    # if np.max(tmp) - np.min(tmp) < 0.02:
    #     Q = R - 1
    from_R = R - Q  # Relative to R
    return rpeak - from_R, Q


def _ecg_delineator_peak_S(rpeak, heartbeat):
    """
    定位S波位置

    Parameters
    ----------
    rpeak
    heartbeat

    Returns
    -------

    """
    segment = heartbeat[rpeak:]  # Select right hand side
    S = signal_findpeaks(-segment, height_min=0.05 * (np.max(segment) - np.min(segment)))

    if len(S["Peaks"]) == 0:
        return math.nan, None
    S = S["Peaks"][0]  # Select most left-hand side
    return rpeak + S, S


def _ecg_delineator_peak_T(rpeak, heartbeat, R, S):
    """
    定位T波起止位置

    Parameters
    ----------
    rpeak
    heartbeat
    R
    S

    Returns
    -------

    """
    if S is None:
        return math.nan, None

    segment = heartbeat[R + S:]  # Select right of S wave
    T = signal_findpeaks(segment, height_min=0.05 * (np.max(segment) - np.min(segment)))

    if len(T["Peaks"]) == 0:
        return math.nan, None
    T = S + T["Peaks"][np.argmax(T["Height"])]  # Select heighest
    return rpeak + T, T


def _ecg_delineator_peak_P_onset(rpeak, heartbeat, R, P):
    """
    P波起点定位

    Parameters
    ----------
    rpeak
    heartbeat
    R
    P

    Returns
    -------

    """
    if P is None:
        return math.nan

    segment = heartbeat.iloc[:P]  # Select left of P wave
    try:
        signal = signal_smooth(segment["Signal"].values, size=R / 10)
    except TypeError:
        signal = segment["Signal"]

    if len(signal) < 2:
        return math.nan

    signal = np.gradient(np.gradient(signal))
    P_onset = np.argmax(signal)

    from_R = R - P_onset  # Relative to R
    return rpeak - from_R


def _ecg_delineator_peak_T_offset(rpeak, heartbeat, R, T):
    """
    T波结束点定位

    Parameters
    ----------
    rpeak
    heartbeat
    R
    T

    Returns
    -------

    """
    if T is None:
        return math.nan

    segment = heartbeat.iloc[R + T:]  # Select right of T wave
    try:
        signal = signal_smooth(segment["Signal"].values, size=R / 10)
    except TypeError:
        signal = segment["Signal"]

    if len(signal) < 2:
        return math.nan

    signal = np.gradient(np.gradient(signal))
    T_offset = np.argmax(signal)

    return rpeak + T + T_offset


def _dwt_compute_multiscales(ecg: np.ndarray, max_degree):
    """Return multiscales wavelet transforms."""

    def _apply_H_filter(signal_i, power=0):
        zeros = np.zeros(2 ** power - 1)
        timedelay = 2 ** power
        banks = np.r_[
            1.0 / 8,
            zeros,
            3.0 / 8,
            zeros,
            3.0 / 8,
            zeros,
            1.0 / 8,
        ]
        signal_f = scipy.signal.convolve(signal_i, banks, mode="full")
        signal_f[:-timedelay] = signal_f[timedelay:]  # timeshift: 2 steps
        return signal_f

    def _apply_G_filter(signal_i, power=0):
        zeros = np.zeros(2 ** power - 1)
        timedelay = 2 ** power
        banks = np.r_[2, zeros, -2]
        signal_f = scipy.signal.convolve(signal_i, banks, mode="full")
        signal_f[:-timedelay] = signal_f[timedelay:]  # timeshift: 1 step
        return signal_f

    dwtmatr = []
    intermediate_ret = np.array(ecg)
    for deg in range(max_degree):
        S_deg = _apply_G_filter(intermediate_ret, power=deg)
        T_deg = _apply_H_filter(intermediate_ret, power=deg)
        dwtmatr.append(S_deg)
        intermediate_ret = np.array(T_deg)
    dwtmatr = [arr[: len(ecg)] for arr in dwtmatr]  # rescale transforms to the same length
    return np.array(dwtmatr)


def dwt_ecg_delineator(heartbeat, rpeak, sampling_rate, analysis_sampling_rate=2000):
    """
    单个Beat的各波定位

    Parameters
    ----------
    heartbeat : 平均波形
    rpeak : 平均波形的R波位置
    sampling_rate : 平均波形采样率
    analysis_sampling_rate : 分析采样率

    Returns
    -------
    P波、T波起始、峰值及结束位置；
    QRS波起止位置
    Q波和S波峰值位置
    """
    # Get index of R peaks
    R = rpeak
    # Q wave
    Q_index, Q = _ecg_delineator_peak_Q(rpeak, heartbeat, R, sampling_rate)
    if math.isnan(Q_index):
        Q_index = int(R - 0.02 * sampling_rate)
    qpeak = Q_index
    resampled_q_index = int(Q_index / sampling_rate * analysis_sampling_rate)
    # S wave
    S_index, S = _ecg_delineator_peak_S(rpeak, heartbeat)
    if math.isnan(S_index):
        S_index = int(R + 0.02 * sampling_rate)
    speak = S_index
    resampled_s_index = int(S_index / sampling_rate * analysis_sampling_rate)

    # dwt to delineate tp waves, onsets, offsets and qrs ontsets and offsets
    heartbeat = resample_interp(heartbeat, sampling_rate, analysis_sampling_rate)
    dwtmatr = _dwt_compute_multiscales(heartbeat, 9)

    rpeak_resampled = int(rpeak / sampling_rate * analysis_sampling_rate)
    average_rate = int(round(60 / (len(heartbeat) / analysis_sampling_rate)))
    tpeak, ppeak = _dwt_delineate_tp_peaks(heartbeat, rpeak_resampled, dwtmatr, average_rate,
                                           sampling_rate=analysis_sampling_rate)
    qrs_onset, qrs_offset = _dwt_delineate_qrs_bounds(rpeak_resampled, dwtmatr, ppeak, tpeak, average_rate,
                                                      resampled_q_index, resampled_s_index,
                                                      sampling_rate=analysis_sampling_rate)

    if math.isnan(qrs_onset):
        if not math.isnan(Q_index):
            qrs_onset = int(resampled_q_index - 0.05 * analysis_sampling_rate)
        else:
            if not math.isnan(ppeak):
                qrs_onset = int(rpeak_resampled - 0.5 * (rpeak_resampled - ppeak))
            else:
                qrs_onset = int(rpeak_resampled - 0.5 * analysis_sampling_rate)
    if math.isnan(qrs_offset):
        if not math.isnan(S_index):
            qrs_offset = int(resampled_s_index - 0.05 * analysis_sampling_rate)
        else:
            if not math.isnan(tpeak):
                qrs_offset = int(rpeak_resampled + 0.2 * (rpeak_resampled - ppeak))
            else:
                qrs_offset = int(rpeak_resampled + 0.5 * analysis_sampling_rate)
    # P bounds
    ponset, poffset = _dwt_delineate_tp_onsets_offsets(ppeak, dwtmatr, average_rate, [qrs_onset, qrs_offset],
                                                       sampling_rate=analysis_sampling_rate, degree_onset=2,
                                                       degree_offset=-3)
    if not math.isnan(ppeak):
        if math.isnan(ponset):
            ponset = max(int(ppeak - 0.15 * analysis_sampling_rate), 0)
        if math.isnan(poffset):
            poffset = min(int(ppeak + 0.15 * analysis_sampling_rate), qrs_onset)
    if qrs_onset - poffset < 0.04 * analysis_sampling_rate:
        poffset = int(ppeak + 0.4 * (qrs_onset - ppeak))

    # T bounds
    tonset, toffset = _dwt_delineate_tp_onsets_offsets(tpeak, dwtmatr, average_rate, [qrs_onset, qrs_offset],
                                                       sampling_rate=analysis_sampling_rate, onset_weight=0.6,
                                                       duration_onset=0.6, degree_onset=2,
                                                       degree_offset=2)
    if not math.isnan(tpeak):
        if math.isnan(tonset):
            tonset = max(int(tpeak - 0.2 * analysis_sampling_rate), 0)
        if math.isnan(toffset):
            toffset = min(int(tpeak + 0.2 * analysis_sampling_rate), int(len(heartbeat) * 0.9))

    if not math.isnan(tpeak) and not math.isnan(ppeak):
        if tonset - qrs_offset < 0.04 * analysis_sampling_rate:
            tonset = int(ppeak + 0.3 * (tpeak - qrs_offset))

    return dict(
        ECG_P_Peak=ppeak if math.isnan(ppeak) else int(ppeak / analysis_sampling_rate * sampling_rate),
        ECG_P_Onset=ponset if math.isnan(ponset) else int(ponset / analysis_sampling_rate * sampling_rate),
        ECG_P_Offset=poffset if math.isnan(poffset) else int(poffset / analysis_sampling_rate * sampling_rate),
        ECG_Q_Peak=qpeak if math.isnan(qpeak) else int(qpeak),
        ECG_R_Onset=qrs_onset if math.isnan(qrs_onset) else int(qrs_onset / analysis_sampling_rate * sampling_rate),
        ECG_R_Offset=qrs_offset if math.isnan(qrs_offset) else int(qrs_offset / analysis_sampling_rate * sampling_rate),
        ECG_S_Peak=speak if math.isnan(speak) else int(speak),
        ECG_T_Peak=tpeak if math.isnan(tpeak) else int(tpeak / analysis_sampling_rate * sampling_rate),
        ECG_T_Onset=tonset if math.isnan(tonset) else int(tonset / analysis_sampling_rate * sampling_rate),
        ECG_T_Offset=toffset if math.isnan(toffset) else int(toffset / analysis_sampling_rate * sampling_rate),
    )
