"""
提取QRS电轴数据
"""
import math


def calculate_axis_using_i_iii_lead(q_r_s_amplitude_i, q_r_s_amplitude_iii):
    """
    通过I和III导联的Q、R、S的振幅和，计算QRS电轴

    Parameters
    ----------
    q_r_s_amplitude_i : I导联上Q、R、S的振幅和
    q_r_s_amplitude_iii : III导联上Q、R、S的振幅和

    Returns
    -------
    电轴度数及偏转类型

    """
    # 电轴
    axis = None
    deviation_type = None
    # 估计avF导联电轴
    vector_avf = q_r_s_amplitude_iii + q_r_s_amplitude_i / 2
    # 排除±90度
    if q_r_s_amplitude_i == 0:
        if vector_avf > 0:
            axis = 90
            deviation_type = 'Normal'
        elif vector_avf < 0:
            axis = -90
            deviation_type = 'Left'
        return axis, deviation_type

    # 计算度数
    x = abs(2 / math.sqrt(3) * (vector_avf / q_r_s_amplitude_i))
    axis_tmp = math.degrees(math.atan(x))
    # 调整度数分区
    if q_r_s_amplitude_i > 0 and vector_avf > 0:
        axis = int(round(axis_tmp))
    elif q_r_s_amplitude_i > 0 and vector_avf < 0:
        axis = -int(round(axis_tmp))
    elif q_r_s_amplitude_i < 0 and vector_avf > 0:
        axis = 180 - int(round(axis_tmp))
    elif q_r_s_amplitude_i < 0 and vector_avf < 0:
        axis = -180 + int(round(axis_tmp))
    if axis is not None:
        if -30 <= axis <= 90:
            deviation_type = 'Normal'
        elif 90 < axis <= 180:
            deviation_type = 'Right'
        elif -90 <= axis < -30:
            deviation_type = 'Left'
        elif -180 < axis < -90:
            deviation_type = 'Extreme'
    return axis, deviation_type


def calculate_q_r_s_amplitude(q_amplitude, r_amplitude, s_amplitude):
    """
    计算综合振幅

    Parameters
    ----------
    q_amplitude : Q波振幅
    r_amplitude : R波振幅
    s_amplitude : S波振幅

    Returns
    -------
    综合振幅

    """
    if q_amplitude is None:
        q_amplitude = 0

    if r_amplitude is None:
        r_amplitude = 0

    if s_amplitude is None:
        s_amplitude = 0

    return r_amplitude - q_amplitude - s_amplitude
