import math

import numpy as np
import scipy.signal as sps
from scipy.signal import find_peaks
from tools.SimpleFilter import smooth_avg1


def _one_peak(dat, index, base_line, qrs_end, qrs_start, sample_rate, qrs):
    if qrs_start < index < qrs_end:
        amplitude = dat[index] - base_line
        if amplitude > 0:
            qrs['r1_index'] = int(index)
            qrs['r1_amplitude'] = float(abs(amplitude))
            qrs['type'] = 'R'
            qrs['r1_duration'] = int((qrs_end - qrs_start) / sample_rate * 1000)
        elif amplitude < 0:
            # Q点
            qrs['q_index'] = int(index)
            qrs['q_amplitude'] = float(abs(amplitude))
            qrs['s_index'] = int(index)
            qrs['s_amplitude'] = float(abs(amplitude))
            qrs['type'] = 'QS'
            qrs['q_duration'] = int((qrs_end - index) / sample_rate * 1000)
            qrs['s_duration'] = int((index - qrs_start) / sample_rate * 1000)
    else:
        Warning('illegal index {}'.format(index))


def _two_peaks(dat, index1, index2, base_line, qrs_end, qrs_start, sample_rate, qrs):
    if index1 != index2:
        tmp = max(index1, index2)
        index1 = min(index1, index2)
        index2 = tmp
        index1_amplitude = dat[index1] - base_line
        index2_amplitude = dat[index2] - base_line
        if index1_amplitude > 0 and index2_amplitude > 0:
            # 转换为R类型
            idx = index1 if abs(index1_amplitude) >= abs(index2_amplitude) else index2
            _one_peak(dat, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)
        elif index1_amplitude < 0 and index2_amplitude < 0:
            # 转换为QS型
            idx = index1 if abs(index1_amplitude) >= abs(index2_amplitude) else index2
            _one_peak(dat, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)
        elif index1_amplitude == 0 or index2_amplitude == 0:
            # 找到的至少一个振幅为0
            # 振幅全为0
            if index1_amplitude == 0 and index2_amplitude == 0:
                Warning('peaks\' amplitude is zero.{}'.format([index1, index2]))
            else:
                # 有一个不为0
                idx = 0 if index1_amplitude == 0 else 1
                index = index2 if idx else index1
                _one_peak(dat, index, base_line, qrs_end, qrs_start, sample_rate, qrs)
        elif index1_amplitude > 0 > index2_amplitude:
            # +-
            # 含有R和S波
            qrs['r1_amplitude'] = float(abs(index1_amplitude))
            qrs['r1_index'] = int(index1)
            qrs['s_amplitude'] = float(abs(index2_amplitude))
            qrs['s_index'] = int(index2)
            total = float(abs(index1_amplitude) + abs(index2_amplitude))
            diff = index2 - index1
            qrs['r1_duration'] = int(
                (abs(index1_amplitude) / total * diff + index1 - qrs_start) / sample_rate * 1000)
            qrs['s_duration'] = int((abs(index2_amplitude) / total * diff + qrs_end - index2) / sample_rate * 1000)
            max_ = max(abs(index1_amplitude), abs(index2_amplitude))
            type_ = ''
            if abs(index1_amplitude) >= max_:
                type_ += 'R'
            else:
                type_ += 'r'

            if abs(index2_amplitude) >= max_:
                type_ += 'S'
            else:
                type_ += 's'
            qrs['type'] = type_
        else:
            # -+
            # 含有Q和R波
            qrs['q_amplitude'] = float(abs(index1_amplitude))
            qrs['q_index'] = int(index1)
            qrs['r1_amplitude'] = float(abs(index2_amplitude))
            qrs['r1_index'] = int(index2)
            total = abs(index1_amplitude) + abs(index2_amplitude)
            diff = int(index2 - index1)
            qrs['q_duration'] = int(
                (abs(index1_amplitude) / total * diff + index1 - qrs_start) / sample_rate * 1000)
            qrs['r1_duration'] = int(
                (abs(index2_amplitude) / total * diff + qrs_end - index2) / sample_rate * 1000)
            max_ = max(abs(index1_amplitude), abs(index2_amplitude))
            type_ = ''
            if abs(index1_amplitude) >= max_:
                type_ += 'Q'
            else:
                type_ += 'q'

            if abs(index2_amplitude) >= max_:
                type_ += 'R'
            else:
                type_ += 'r'
            qrs['type'] = type_
    else:
        _one_peak(dat, index1, base_line, qrs_end, qrs_start, sample_rate, qrs)


def _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, index1, index2, index3):
    indexes = [index1, index2, index3]
    indexes.sort()
    index1 = indexes[0]
    index2 = indexes[1]
    index3 = indexes[2]
    index1_amplitude = dat[index1] - base_line
    index2_amplitude = dat[index2] - base_line
    index3_amplitude = dat[index3] - base_line
    arr = np.array([index1_amplitude, index2_amplitude, index3_amplitude])
    if index1_amplitude == 0 or index2_amplitude == 0 or index3_amplitude == 0:
        not_zero_indexes = np.array(indexes)[np.array(arr != 0)]
        if len(not_zero_indexes) == 0:
            Warning('empty peaks')
        elif len(not_zero_indexes) == 1:
            _one_peak(dat, not_zero_indexes[0], base_line, qrs_end, qrs_start, sample_rate, qrs)
        else:
            _two_peaks(dat, not_zero_indexes[0], not_zero_indexes[1], base_line, qrs_end, qrs_start, sample_rate, qrs)
    else:
        # 3个点的振幅没有为0的
        if index1_amplitude > 0:
            # +xx
            # 第一个点振幅大于0
            if index2_amplitude > 0:
                # ++x
                if index3_amplitude > 0:
                    # +++
                    max_ = np.max(arr)
                    idx = np.array([index1, index2, index3])[np.where(arr == max_)[0][0]]
                    _one_peak(dat, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)
                else:
                    # ++-
                    (idx_1, idx_2) = (index1, index3) if index1_amplitude >= index2_amplitude else (index2, index3)
                    _two_peaks(dat, idx_1, idx_2, base_line, qrs_end, qrs_start, sample_rate, qrs)
            else:
                # +-x
                if index3_amplitude > 0:
                    # +-+
                    # rsr'类型
                    qrs['r1_index'] = index1
                    qrs['s_index'] = index2
                    qrs['r2_index'] = index3
                    qrs['r1_amplitude'] = abs(index1_amplitude)
                    qrs['s_amplitude'] = abs(index2_amplitude)
                    qrs['r2_amplitude'] = abs(index3_amplitude)
                    total_1 = abs(index2_amplitude) + abs(index1_amplitude)
                    diff_1 = index2 - index1
                    total_2 = abs(index2_amplitude) + abs(index3_amplitude)
                    diff_2 = index3 - index2
                    qrs['r1_duration'] = int(
                        (index1 - qrs_start + qrs['r1_amplitude'] / total_1 * diff_1) / sample_rate * 1000)
                    qrs['s_duration'] = int(
                        (qrs['s_amplitude'] / total_2 * diff_2 + qrs[
                            's_amplitude'] / total_1 * diff_1) / sample_rate * 1000)
                    qrs['r2_duration'] = int(
                        (qrs_end - index2 + qrs['r2_amplitude'] / total_2 * diff_2) / sample_rate * 1000)
                    type_ = ''
                    max_ = np.max(np.abs(arr))
                    if abs(index1_amplitude) >= max_:
                        type_ += 'R'
                    else:
                        type_ += 'r'

                    if abs(index2_amplitude) >= max_:
                        type_ += 'S'
                    else:
                        type_ += 's'

                    if abs(index3_amplitude) >= max_:
                        type_ += 'R\''
                    else:
                        type_ += 'r\''

                    qrs['type'] = type_
                else:
                    # +--
                    (idx_1, idx_2) = (index1, index2) if index2_amplitude < index2_amplitude else (index1, index3)
                    _two_peaks(dat, idx_1, idx_2, base_line, qrs_end, qrs_start, sample_rate, qrs)
        else:
            # -xx
            # 第一个点振幅小于0
            if index2_amplitude > 0:
                # -+x
                if index3_amplitude > 0:
                    # -++
                    (idx_1, idx_2) = (index1, index2) if index2_amplitude >= index3_amplitude else (index1, index3)
                    _two_peaks(dat, idx_1, idx_2, base_line, qrs_end, qrs_start, sample_rate, qrs)
                else:
                    # -+-
                    # qrs类型
                    qrs['q_index'] = index1
                    qrs['r1_index'] = index2
                    qrs['s_index'] = index3
                    qrs['q_amplitude'] = float(abs(index1_amplitude))
                    qrs['r1_amplitude'] = float(abs(index2_amplitude))
                    qrs['s_amplitude'] = float(abs(index3_amplitude))
                    total_1 = abs(index2_amplitude) + abs(index1_amplitude)
                    diff_1 = index2 - index1
                    total_2 = abs(index2_amplitude) + abs(index3_amplitude)
                    diff_2 = index3 - index2
                    qrs['q_duration'] = int(
                        (index1 - qrs_start + qrs['q_amplitude'] / total_1 * diff_1) / sample_rate * 1000)
                    qrs['r1_duration'] = int(
                        (qrs['r1_amplitude'] / total_2 * diff_2 + qrs[
                            'r1_amplitude'] / total_1 * diff_1) / sample_rate * 1000)
                    qrs['s_duration'] = int(
                        (qrs_end - index2 + qrs['s_amplitude'] / total_2 * diff_2) / sample_rate * 1000)
                    type_ = ''
                    max_ = np.max(np.abs(arr))
                    if abs(index1_amplitude) >= max_:
                        type_ += 'Q'
                    else:
                        type_ += 'q'

                    if abs(index2_amplitude) >= max_:
                        type_ += 'R'
                    else:
                        type_ += 'r'

                    if abs(index3_amplitude) >= max_:
                        type_ += 'S'
                    else:
                        type_ += 's'

                    qrs['type'] = type_

            else:
                # --x
                if index3_amplitude > 0:
                    # --+
                    (idx_1, idx_2) = (index2, index3) if index1_amplitude > index2_amplitude else (index1, index3)
                    _two_peaks(dat, idx_1, idx_2, base_line, qrs_end, qrs_start, sample_rate, qrs)
                else:
                    # ---
                    min_ = np.min(arr)
                    idx = np.array([index1, index2, index3])[np.where(arr == min_)[0][0]]
                    _one_peak(dat, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)


def _four_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, index1, index2, index3, index4):
    indexes = [index1, index2, index3, index4]
    indexes.sort()
    index1 = indexes[0]
    index2 = indexes[1]
    index3 = indexes[2]
    index4 = indexes[3]
    index1_amplitude = dat[index1] - base_line
    index2_amplitude = dat[index2] - base_line
    index3_amplitude = dat[index3] - base_line
    index4_amplitude = dat[index4] - base_line
    arr = np.array([index1_amplitude, index2_amplitude, index3_amplitude, index4_amplitude])
    if len(np.where(arr == 0)[0]) > 0:
        not_zero_indexes = np.array(indexes)[np.array(arr != 0)]
        if len(not_zero_indexes) == 0:
            Warning('empty peaks')
        elif len(not_zero_indexes) == 1:
            _one_peak(dat, not_zero_indexes[0], base_line, qrs_end, qrs_start, sample_rate, qrs)
        elif len(not_zero_indexes) == 2:
            _two_peaks(dat, not_zero_indexes[0], not_zero_indexes[1], base_line, qrs_end, qrs_start, sample_rate, qrs)
        elif len(not_zero_indexes) == 3:
            _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, not_zero_indexes[0], not_zero_indexes[1],
                         not_zero_indexes[2])
    else:
        if index1_amplitude > 0:
            # +xxx
            if index2_amplitude > 0:
                # ++xx
                if index3_amplitude > 0:
                    # +++x
                    idx = np.array([index1, index2, index3, index4])[np.where(arr == np.max(arr))[0][0]]
                    if index4_amplitude > 0:
                        # ++++
                        _one_peak(dat, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)
                    else:
                        # +++-
                        _two_peaks(dat, idx, index4, base_line, qrs_end, qrs_start, sample_rate, qrs)
                else:
                    # ++-x
                    if index4_amplitude > 0:
                        # ++-+
                        idx = index1 if index1_amplitude >= index2_amplitude else index2
                        _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, idx, index3, index4)
                    else:
                        # ++--
                        idx1 = index1 if index1_amplitude >= index2_amplitude else index2
                        idx2 = index3 if abs(index3_amplitude) >= abs(index4_amplitude) else index4
                        _two_peaks(dat, idx1, idx2, base_line, qrs_end, qrs_start, sample_rate, qrs)
            else:
                # +-xx
                if index3_amplitude > 0:
                    # +-+x
                    if index4_amplitude > 0:
                        # +-++
                        idx = index3 if index3_amplitude >= index4_amplitude else index4
                        _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, index1, index2, idx)
                    else:
                        # +-+-
                        qrs['r1_index'] = index1
                        qrs['s_index'] = index2
                        qrs['r2_index'] = index3
                        qrs['s2_index'] = index4
                        qrs['r1_amplitude'] = float(abs(index1_amplitude))
                        qrs['s_amplitude'] = float(abs(index2_amplitude))
                        qrs['r2_amplitude'] = float(abs(index3_amplitude))
                        qrs['s2_amplitude'] = float(abs(index4_amplitude))
                        total1 = abs(index1_amplitude) + abs(index2_amplitude)
                        diff1 = index2 - index1
                        total2 = abs(index2_amplitude) + abs(index3_amplitude)
                        diff2 = index3 - index2
                        total3 = abs(index3_amplitude) + abs(index4_amplitude)
                        diff3 = index4 - index3
                        qrs['r1_duration'] = int(
                            (index1 - qrs_start + abs(index1_amplitude) / total1 * diff1) / sample_rate * 1000)
                        qrs['s_duration'] = int(
                            (abs(index2_amplitude) / total2 * diff2 + abs(
                                index2_amplitude) / total1 * diff1) / sample_rate * 1000)
                        qrs['r2_duration'] = int(
                            (abs(index3_amplitude) / total2 * diff2 + abs(
                                index3_amplitude) / total3 * diff3) / sample_rate * 1000)
                        qrs['s2_duration'] = int(
                            (abs(index4_amplitude) / total3 * diff3 + qrs_end - index4) / sample_rate * 1000)

                        type_ = ''
                        max_ = np.max(np.abs(arr))
                        if abs(index1_amplitude) >= max_:
                            type_ += 'R'
                        else:
                            type_ += 'r'

                        if abs(index2_amplitude) >= max_:
                            type_ += 'S'
                        else:
                            type_ += 's'

                        if abs(index3_amplitude) >= max_:
                            type_ += 'R\''
                        else:
                            type_ += 'r\''

                        if abs(index4_amplitude) >= max_:
                            type_ += 'S\''
                        else:
                            type_ += 's\''

                        qrs['type'] = type_

                else:
                    # +--x
                    if index4_amplitude > 0:
                        # +--+
                        idx = index2 if abs(index2) >= abs(index3) else index3
                        _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, index1, idx, index4)
                    else:
                        # +---
                        idx = np.array([index1, index2, index3, index4])[np.where(arr == np.min(arr))[0][0]]
                        _two_peaks(dat, index1, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)
        else:
            # -xxx
            if index2_amplitude > 0:
                # -+xx
                if index3_amplitude > 0:
                    # -++x
                    idx = np.array([index1, index2, index3, index4])[np.where(arr == np.max(arr))[0][0]]
                    if index4_amplitude > 0:
                        # -+++
                        _two_peaks(dat, index1, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)
                    else:
                        # -++-
                        _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, index1, idx, index4)
                else:
                    # -+-x
                    if index4_amplitude > 0:
                        # -+-+
                        qrs['q_index'] = int(index1)
                        qrs['r1_index'] = int(index2)
                        qrs['s_index'] = int(index3)
                        qrs['r2_index'] = int(index4)
                        qrs['q_amplitude'] = float(abs(index1_amplitude))
                        qrs['r1_amplitude'] = float(abs(index2_amplitude))
                        qrs['s_amplitude'] = float(abs(index3_amplitude))
                        qrs['r2_amplitude'] = float(abs(index4_amplitude))
                        total1 = float(abs(index1_amplitude) + abs(index2_amplitude))
                        diff1 = index2 - index1
                        total2 = float(abs(index2_amplitude) + abs(index3_amplitude))
                        diff2 = index3 - index2
                        total3 = float(abs(index3_amplitude) + abs(index4_amplitude))
                        diff3 = index4 - index3
                        qrs['q_duration'] = int(
                            (index1 - qrs_start + abs(index1_amplitude) / total1 * diff1) / sample_rate * 1000)
                        qrs['r1_duration'] = int(
                            (abs(index2_amplitude) / total2 * diff2 + abs(
                                index2_amplitude) / total1 * diff1) / sample_rate * 1000)
                        qrs['s_duration'] = int(
                            (abs(index3_amplitude) / total2 * diff2 + abs(
                                index3_amplitude) / total3 * diff3) / sample_rate * 1000)
                        qrs['r2_duration'] = int(
                            (abs(index4_amplitude) / total3 * diff3 + qrs_end - index4) / sample_rate * 1000)

                        type_ = ''
                        max_ = np.max(np.abs(arr))
                        if abs(index1_amplitude) >= max_:
                            type_ += 'Q'
                        else:
                            type_ += 'q'

                        if abs(index2_amplitude) >= max_:
                            type_ += 'R'
                        else:
                            type_ += 'r'

                        if abs(index3_amplitude) >= max_:
                            type_ += 'S'
                        else:
                            type_ += 's'

                        if abs(index4_amplitude) >= max_:
                            type_ += 'R\''
                        else:
                            type_ += 'r\''

                        qrs['type'] = type_
                    else:
                        # -+--
                        idx = index3 if abs(index3_amplitude) >= abs(index4_amplitude) else index4
                        _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, index1, index2, idx)
            else:
                # --xx
                if index3_amplitude > 0:
                    # --+x
                    if index4_amplitude > 0:
                        # --++
                        idx1 = index1 if abs(index1_amplitude) >= abs(index2_amplitude) else index2
                        idx2 = index3 if abs(index3_amplitude) >= abs(index4_amplitude) else index4
                        _two_peaks(dat, idx1, idx2, base_line, qrs_end, qrs_start, sample_rate, qrs)
                    else:
                        # --+-
                        idx1 = index1 if abs(index1_amplitude) >= abs(index2_amplitude) else index2
                        _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, idx1, index3, index4)
                else:
                    # ---x
                    idx = np.array([index1, index2, index3, index4])[np.where(arr == np.min(arr))[0][0]]
                    if index4_amplitude > 0:
                        # ---+
                        _two_peaks(dat, idx, index4, base_line, qrs_end, qrs_start, sample_rate, qrs)
                    else:
                        # ----
                        _one_peak(dat, idx, base_line, qrs_end, qrs_start, sample_rate, qrs)


def _all_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, indexes):
    indexes = list(set(indexes))
    indexes.sort()
    if len(indexes) == 4:
        _four_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, indexes[0],
                    indexes[1], indexes[2], indexes[3])
    elif len(indexes) == 3:
        _three_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, indexes[0],
                     indexes[1], indexes[2])
    elif len(indexes) == 2:
        _two_peaks(dat, indexes[0], indexes[1], base_line, qrs_end, qrs_start, sample_rate, qrs)
    elif len(indexes) == 1:
        _one_peak(dat, indexes[0], base_line, qrs_end, qrs_start, sample_rate, qrs)
    else:
        amplitudes = []
        for _ in indexes:
            amplitudes.append(dat[_] - base_line)
        amplitudes = np.array(amplitudes)
        amplitudes_sort = np.sort(np.array(amplitudes))
        indexes_ = []
        max_index = np.where(amplitudes == amplitudes_sort[-1])[0]
        if len(max_index) > 1:
            indexes_.append(max_index[-1])
            indexes_.append(max_index[-2])
        else:
            indexes_.append(max_index[-1])
            max_index = np.where(amplitudes == amplitudes_sort[-2])[0]
            indexes_.append(max_index[-1])

        min_index = np.where(amplitudes == amplitudes_sort[0])[0]
        if len(min_index) > 1:
            indexes_.append(min_index[-1])
            indexes_.append(min_index[-2])
        else:
            indexes_.append(min_index[-1])
            min_index = np.where(amplitudes == amplitudes_sort[1])[0]
            indexes_.append(min_index[-1])
        _four_peaks(dat, qrs, qrs_end, qrs_start, base_line, sample_rate, indexes[indexes_[0]],
                    indexes[indexes_[1]], indexes[indexes_[2]], indexes[indexes_[3]])


def _qrs_complex(baseline, qrs_end, qrs_start, sample_rate, wave_data):
    qrs = {'duration': int((qrs_end - qrs_start) / sample_rate * 1000)}
    offset = int(0.1 * sample_rate)
    qrs_wave_data = wave_data[max(int(qrs_start - offset), 0):int(qrs_end + offset)]
    # 峰值检测
    smooth_qrs_dat = smooth_avg1(qrs_wave_data, 4)
    peaks = find_peaks(smooth_qrs_dat)[0] + qrs_start - offset
    peaks = peaks[np.where(peaks > qrs_start)]
    peaks = peaks[np.where(peaks < qrs_end)]
    result_peaks = []
    for _ in peaks:
        if abs(_ - qrs_start) <= 0.02 * sample_rate:
            if wave_data[_] - wave_data[qrs_start] > 0.05:
                result_peaks.append(_)
        elif abs(_ - qrs_end) <= 0.02 * sample_rate:
            if wave_data[_] - wave_data[qrs_end] > 0.05:
                result_peaks.append(_)
        else:
            result_peaks.append(_)

    peaks = find_peaks(-np.array(smooth_qrs_dat))[0] + qrs_start - offset
    peaks = peaks[np.where(peaks > qrs_start)]
    peaks = peaks[np.where(peaks < qrs_end)]
    for _ in peaks:
        if abs(_ - qrs_start) <= 0.02 * sample_rate:
            if wave_data[_] - wave_data[qrs_start] > 0.05:
                result_peaks.append(_)
        elif abs(_ - qrs_end) <= 0.02 * sample_rate:
            if wave_data[_] - wave_data[qrs_end] > 0.05:
                result_peaks.append(_)
        else:
            result_peaks.append(_)
    if len(result_peaks) < 1:
        result_peaks.append(np.argmax(smooth_qrs_dat[1:-2]))

    _all_peaks(wave_data, qrs, qrs_end, qrs_start, baseline, sample_rate, result_peaks)
    r1 = qrs.get('r1_amplitude', 0)
    r2 = qrs.get('r2_amplitude', 0)
    q = qrs.get('q_amplitude', 0)
    s = qrs.get('s_amplitude', 0)
    if max(r1, r2) >= max(q, s):
        qrs['qrs_principle'] = '+'
    else:
        qrs['qrs_principle'] = '-'
    qrs['qrs_amplitude'] = max(r1, r2) + max(q, s)

    if q != 0 and (r1 != 0 or r2 != 0):
        qrs['Q/R'] = q / max(r1, r2)
    else:
        qrs['Q/R'] = math.nan

    J = False
    if qrs.get('type', None) is not None:
        if str(qrs['type'])[-1] in ['R', 'r']:
            J_data = wave_data[qrs['r1_index']:qrs_end]
            peaks = find_peaks(J_data)[0]
            if peaks is None or len(peaks) < 1:
                pass
            else:
                J_idx = np.where(np.logical_and(qrs['r1_index'] + 2 < peaks, peaks < qrs_end - 2))[0]
                if len(J_idx) < 1:
                    pass
                else:
                    J = True
    qrs['ExistJ'] = J

    key = 'q_amplitude'
    q_a = qrs.get(key, None)
    if q_a is None:
        qrs[key] = math.nan

    key = 'q_duration'
    q_d = qrs.get(key, None)
    if q_d is None:
        qrs[key] = math.nan

    key = 'r1_amplitude'
    q_d = qrs.get(key, None)
    if q_d is None:
        qrs[key] = math.nan

    key = 'r1_duration'
    q_d = qrs.get(key, None)
    if q_d is None:
        qrs[key] = math.nan

    key = 'r2_amplitude'
    q_d = qrs.get(key, None)
    if q_d is None:
        qrs[key] = math.nan

    key = 'r2_duration'
    q_d = qrs.get(key, None)
    if q_d is None:
        qrs[key] = math.nan

    key = 's_amplitude'
    q_d = qrs.get(key, None)
    if q_d is None:
        qrs[key] = math.nan

    key = 's_duration'
    q_d = qrs.get(key, None)
    if q_d is None:
        qrs[key] = math.nan
    for k, v in qrs.items():
        if isinstance(v, dict):
            pass
        elif isinstance(v, np.float):
            qrs[k] = float(v)
        elif isinstance(v, np.int):
            qrs[k] = int(v)

    return qrs


def _p_or_t_wave(baseline, end, start, sample_rate, wave_data: list):
    result = {'exist': True}
    if end is math.nan or start is math.nan:
        result['exist'] = False
        return result
    wave_data = wave_data[start:end]
    # smooth
    smooth_dat = smooth_avg1(wave_data, int(0.008 * sample_rate))
    # peaks
    max_index = np.argmax(smooth_dat)
    min_index = np.argmin(smooth_dat)
    max_ = wave_data[max_index]
    min_ = wave_data[min_index]

    # form
    if max_ >= baseline and min_ >= baseline:
        result['form'] = '+'
        amplitude = [max(abs(max_ - baseline), abs(min_ - baseline))]
        result['amplitude'] = float(sum(amplitude))
        result['direction'] = '+'
    elif max_ <= baseline and min_ <= baseline:
        result['form'] = '-'
        amplitude = [max(abs(max_ - baseline), abs(min_ - baseline))]
        result['amplitude'] = float(sum(amplitude))
        result['direction'] = '-'
    else:
        if max_index > min_index:
            result['form'] = '-+'
            amplitude = [abs(min_ - baseline), abs(max_ - baseline)]
            result['amplitude'] = float(sum(amplitude))
            if amplitude[0] >= amplitude[1]:
                result['direction'] = '-'
                if abs(amplitude[0] - baseline) < 0.05:
                    result['form'] = '-'
            else:
                result['direction'] = '+'
                if abs(amplitude[1] - baseline) < 0.05:
                    result['form'] = '+'
        else:
            result['form'] = '+-'
            amplitude = [abs(max_ - baseline), abs(min_ - baseline)]
            result['amplitude'] = float(sum(amplitude))
            if amplitude[0] >= amplitude[1]:
                result['direction'] = '+'
                if abs(amplitude[1] - baseline) < 0.05:
                    result['form'] = '+'
            else:
                result['direction'] = '-'
                if abs(amplitude[0] - baseline) < 0.05:
                    result['form'] = '-'
    result['duration'] = int(round((end - start) * 1000 / sample_rate))
    return result


def _pr_or_st(wave_data, fs, start1, end1, baseline, start2=None):
    res = {'computable': True}
    if math.isnan(start1) or math.isnan(end1) or end1 < start1:
        res['computable'] = False
        return res
    if start2 is not None and math.isnan(start2):
        res['computable'] = False
        return res
    if start2 is not None:
        res['interval'] = int((end1 - start2) / fs * 1000)
    else:
        res['duration'] = int((end1 - start1) / fs * 1000)
    # PR form
    tmp_data = wave_data[start1-2:end1+2]
    max_v = np.max(tmp_data)
    max_idx = np.argmax(tmp_data)
    min_v = np.min(tmp_data)
    min_idx = np.argmin(tmp_data)
    if abs(max_v - min_v) < 0.05:
        res['form'] = 'horizontal'
        if start2 is None:
            res['amplitude'] = 0
    else:
        if 1 < max_idx < len(tmp_data) - 2 and 1 < min_idx < len(tmp_data) - 2:
            idx = max_idx if abs(max_v - baseline) >= abs(min_v - baseline) else min_idx
            if tmp_data[0] - tmp_data[idx] >= 0:
                res['form'] = 'declination'
            else:
                res['form'] = 'upslope'
            if start2 is None:
                res['amplitude'] = float(tmp_data[idx] - tmp_data[0])

        elif 1 < max_idx < len(tmp_data) - 2:
            idx = max_idx
            if tmp_data[0] - tmp_data[idx] >= 0:
                res['form'] = 'declination'
            else:
                res['form'] = 'upslope'

            if start2 is None:
                res['amplitude'] = float(tmp_data[idx] - tmp_data[0])

        elif 1 < min_idx < len(tmp_data) - 2:
            idx = min_idx
            if tmp_data[0] - tmp_data[idx] >= 0:
                res['form'] = 'declination'
            else:
                res['form'] = 'upslope'

            if start2 is None:
                res['amplitude'] = float(tmp_data[idx] - tmp_data[0])

        else:
            idx = -1
            if tmp_data[0] - tmp_data[idx] >= 0:
                res['form'] = 'declination'
            else:
                res['form'] = 'upslope'

            if start2 is None:
                res['amplitude'] = float(tmp_data[idx] - tmp_data[0])
    return res


def extract_features(wave_data, p_start, p_end, qrs_start, qrs_end, t_start, t_end, sample_rate, avg_rr,
                     adc_gain=1, adc_zero=0):
    """
    波形参数解析
    波形形态、幅值、时限等
    :param wave_data: 波形数据
    :param p_start: p波开始下标
    :param p_end: p波结束下标
    :param qrs_start: qrs起始下标
    :param qrs_end: qrs结束下标
    :param t_end: t波结束下标
    :param sample_rate: 采样率
    :param adc_gain: 增益
    :param adc_zero: 0点电压
    :return:
    """
    if p_start < 0:
        p_start = 0
    wave_data = np.round((np.array(wave_data) - adc_zero) / adc_gain, 3)
    # 以q波起点作为基线
    baseline = wave_data[qrs_start]
    result = {}
    # P
    p_wave = _p_or_t_wave(baseline, p_end, p_start, sample_rate, wave_data)

    # T
    t_wave = _p_or_t_wave(baseline, t_end, qrs_end, sample_rate, wave_data)

    # PR
    pr_wave = _pr_or_st(wave_data, sample_rate, p_end, qrs_start, baseline, p_start)

    # ST
    st_wave = _pr_or_st(wave_data, sample_rate, qrs_end, t_start, baseline)

    # QRS
    qrs_wave = _qrs_complex(baseline, qrs_end, qrs_start, sample_rate, wave_data)
    if not (math.isnan(t_end) or math.isnan(qrs_end) or qrs_end >= t_end):
        qt = int((t_end - qrs_end) / sample_rate * 1000)
        result['QT'] = qt
        result['QTc'] = int(qt / 1000 / math.sqrt(avg_rr / 1000) * 1000)
    else:
        result['QT'] = math.nan
        result['QTc'] = math.nan
    result['HR'] = int(60 / (avg_rr / 1000))
    r1 = 0 if math.isnan(qrs_wave.get('r1_amplitude', math.nan)) else float(qrs_wave['r1_amplitude'])
    r2 = 0 if math.isnan(qrs_wave.get('r2_amplitude', math.nan)) else float(qrs_wave['r2_amplitude'])
    r = max(r1, r2)
    t = t_wave.get('amplitude', 0)
    if r > 0 and t > 0:
        result['R/T'] = float(r / t)
    else:
        result['R/T'] = math.nan

    _QRS = {'RAmplitude': qrs_wave.get('r1_amplitude', math.nan),
            'SAmplitude': qrs_wave.get('s_amplitude', math.nan),
            'QAmplitude': qrs_wave.get('q_amplitude', math.nan),
            'R\'Amplitude': qrs_wave.get('r2_amplitude', math.nan),
            'RDuration': qrs_wave.get('r1_duration', math.nan),
            'SDuration': qrs_wave.get('s_duration', math.nan),
            'QDuration': qrs_wave.get('q_duration', math.nan),
            'R\'Duration': qrs_wave.get('r2_duration', math.nan),
            'QRSPrinciple': qrs_wave.get('qrs_principle', math.nan),
            'QRSAmplitude': qrs_wave.get('qrs_amplitude', math.nan),
            'QRSDuration': qrs_wave.get('duration', math.nan),
            'Q/R': qrs_wave.get('Q/R', math.nan)}

    result['QRS'] = _QRS
    result['PR'] = pr_wave
    result['ST'] = st_wave
    result['T'] = t_wave
    result['P'] = p_wave
    return result

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
    return np.mean(cache_beats, axis=0), ahead_len, avg_len