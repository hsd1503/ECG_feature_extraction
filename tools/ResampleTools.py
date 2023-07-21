"""
重采样工具
"""
import numpy as np
from scipy.interpolate import interp1d


def resample_interp(ts, fs_in, fs_out):
    """
    基于线性拟合的差值重采样算法
    计算前后点对应的比例进行插值
    :param ts:  单导联数据，一维浮点型数组
    :param fs_in: 原始采样率，整型
    :param fs_out: 目标采样率，整型
    :return: 重采样后的数据
    """
    t = len(ts) / fs_in
    fs_in, fs_out = int(fs_in), int(fs_out)
    if fs_out == fs_in:
        return np.array(ts)
    else:
        x_old = np.linspace(0, 1, num=len(ts), endpoint=True)
        x_new = np.linspace(0, 1, num=int(t * fs_out), endpoint=True)
        y_old = ts
        f = interp1d(x_old, y_old, kind='linear')
        y_new = f(x_new)
        return y_new


def resample_with_amplitude(sig, fs_in, fs_out):
    """
    保留较大振幅的重采样方法
    可以在较低采样率的时候保持一定的高频特性


    Parameters
    ----------
    sig : 原始信号
    fs_in : 原始信号采样率
    fs_out : 重采样后的采样率

    Returns
    -------
    重采样后的数据

    """
    desired_size = int(len(sig) / fs_in * fs_out)
    adap_t, adap_x = adaptive_resampling(sig, desired_size)
    warping_time = np.linspace(0, adap_t[-1], desired_size)
    return warping(adap_x, adap_t, warping_time)


def adaptive_resampling(time_sequence, desired_size, type='original'):
    """
    自适应重采样


    Parameters
    -----
    time_sequence: 待重采样时间序列
    desired_size: 期望重采样后的长度
    type: `original`表示使用原始方式，否者使用新的方式


    Returns
    -----
    重采样后序列时间轴点及对应样本点


    """

    sig = np.array(time_sequence)
    D_Q = np.insert(np.cumsum(np.abs(np.diff(sig))), 0, 0)
    d_Q = D_Q[-1] / desired_size
    if 'original' == type:
        t = [np.where(D_Q >= i * d_Q)[0][0] for i in range(desired_size)]
    else:
        t = []
        i = 0
        idx = 0

        while i < desired_size:
            while idx < len(D_Q) and D_Q[idx] < d_Q * i:
                idx += 1
            t.append(idx)
            i += 1
    t = np.array(t)
    beta = np.array([(D_Q[t[i]] - d_Q * i) / (D_Q[t[i]] - D_Q[t[i] - 1]) for i in range(desired_size)])

    t_hat = t - beta
    x_hat = sig[t] - beta * (sig[t] - sig[t - 1])

    return t_hat, x_hat


def warping(sig, original_time, target_time):
    """
    时间规整，将时间尺度上的点重新规整到期望的时间点上

    Parameters
    -------------------------
    sig: 原始信号采样点
    original_time: 原始信号采样点对应时间
    target_time: 期望重新规整后的时间点


    Returns
    -------------------------
    规整后的采样点

    """
    f = interp1d(original_time, sig, kind='linear')
    return f(target_time)
