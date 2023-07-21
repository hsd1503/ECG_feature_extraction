"""
通用工具
"""
from neurokit2.hrv import hrv_utils
import numpy as np
from neurokit2.signal import signal_interpolate
import math


def hrv_get_rri(peaks: list, sampling_rate=1000, interpolate=False, **kwargs):
    """
    从R波下表中获取rr间隙
    :param peaks: R波下标
    :param sampling_rate: 信号采样率
    :param interpolate: 是否重采样
    :return:
    """
    peaks = hrv_utils._hrv_sanitize_input(peaks)
    if isinstance(peaks, tuple):  # Detect actual sampling rate
        peaks, sampling_rate = peaks[0], peaks[1]

    return _hrv_get_rri(peaks, sampling_rate, interpolate, **kwargs)


def _hrv_get_rri(peaks=None, sampling_rate=1000, interpolate=False, **kwargs):
    rri = np.diff(peaks) / sampling_rate * 1000
    if interpolate is False:
        return rri, sampling_rate

    else:
        # Minimum sampling rate for interpolation
        if sampling_rate < 10:
            sampling_rate = 10

        # Compute length of interpolated heart period signal at requested sampling rate.
        desired_length = int(np.rint(peaks[-1]))

        rri = signal_interpolate(
            peaks[1:],  # Skip first peak since it has no corresponding element in heart_period
            rri,
            x_new=np.arange(desired_length),
            **kwargs
        )
        return rri, sampling_rate


class TimeDomainBasicIndices(object):
    """
    时域基础指标封装类
    """
    # 最大RR间隙
    max_rr = None
    # 最小RR间隙
    min_rr = None
    # 平均RR间隙
    mean_rr = None
    # RR间隙中位数
    median_rr = None
    # 中值绝对偏差
    mad_rr = None
    # 归一化的mad_rr
    mcv_rr = None
    # RR间隙方差
    variance_rr = None
    # RR间隙标准差
    std_rr = None
    # RR间隙波动范围
    amplitude_rr = None
    # RR间隙5%分位数
    percentile_5 = None
    # RR间隙25%分位数
    percentile_25 = None
    # RR间隙75%分位数
    percentile_75 = None
    # RR间隙95%分位数
    percentile_95 = None
    # 75%-25%分位数差值
    percentile_range_75_25 = None
    # 95%-5%分位数差值
    percentile_range_95_5 = None
    # 偏度值
    skew_rr = None
    # 峰度值
    kurtosis_rr = None
    # RR间隙数量
    count_rr = None


class TimeDomainStandardIndices(object):
    # RR间隙标准差，自由度为N-1
    sdnn = None
    # 计算每分钟一段的NN间隙均值的标准差
    sdann_1 = None
    # 计算每分钟片段sdnn值的均值
    sdnni_1 = None
    # 计算每2分钟一段的NN间隙均值的标准差
    sdann_2 = None
    # 计算每2分钟片段sdnn值的均值
    sdnni_2 = None
    # 计算每5分钟一段的NN间隙均值的标准差
    sdann_5 = None
    # 计算每5分钟片段sdnn值的均值
    sdnni_5 = None
    # 连续RR间隙差值的均方根
    rmssd = None
    # 连续RR间隙差值的标准差，自由度为N-1
    sdsd = None
    # 归一化后的sdnn
    cvnn = None
    # 归一化后的rmssd
    cvsd = None
    # 计算连续RR间隙差值中，大于等于50ms的数量百分比
    pNN50 = None
    # 计算连续RR间隙差值中，大于等于20ms的数量百分比
    pNN20 = None
    # RR间隙数量与RR间隙直方图中最大值的比值
    tri_index = None
    # 计算直方图中最小二乘的M和N值的最小值
    tinn = None
    # tinn的m值
    tinn_m = None
    # tinn的n值
    tinn_n = None


class FrequencyDomainIndices(object):
    """
    频域指标封装类
    """

    def __init__(self, ulf_band=(0, 0.0033), vlf_band=(0.0033, 0.04),
                 lf_band=(0.04, 0.15), hf_band=(0.15, 0.4),
                 vhf_band=(0.4, 0.5), psd_method="welch", fs=4):
        """
        初始化
        :param ulf_band: 超低频率带宽
        :param vlf_band: 较低频率带宽
        :param lf_band: 低频率带宽
        :param hf_band: 高频率带宽
        :param vhf_band: 较高频率带宽
        :param psd_method: 功率谱分析方法
        :param fs: RR间隙重采样频率
        """
        self.ulf_band = ulf_band
        self.vlf_band = vlf_band
        self.lf_band = lf_band
        self.hf_band = hf_band
        self.vhf_band = vhf_band
        self.psd_method = psd_method
        self.fs = fs
        # 超低频率功率密度
        self.ulf_power = 0
        # 较低频率功率密度
        self.vlf_power = 0
        # 低频率功率密度
        self.lf_power = 0
        # 高频率功率密度
        self.hf_power = 0
        # 较高频率功率密度
        self.vhf_power = 0
        # 总功率密度
        self.total_power = 0
        # 低频功率密度与总功率密度比值
        self.lf_n = 0
        # 高频功率密度与总功率密度比值
        self.hf_n = 0
        # 低功率密度与高功率密度比值
        self.lf_hf = 0


class NonLinearPoincareIndices(object):
    """
    非线性分析中的庞加莱相关指标
    """
    # 庞加莱图中所有点所在面积
    s = None
    # 庞加莱图中垂直于特征线的长度
    sd1 = None
    # 庞加莱图中特征线的长度
    sd2 = None
    # sd_1与sd_2比值
    sd1_sd2 = None
    # 心脏交感指数
    csi = None
    # 心脏迷走神经指数
    cvi = None
    # 修正后的心脏交感指数
    csi_modified = None


class NonLinearEntropyIndices(object):
    """
    非线性分析中的熵值相关指标
    """
    # 近似熵
    approximate = None

    # 样本熵
    sample = None

    # 香侬熵
    shannon = None

    # 模糊熵
    fuzzy = None

    # 标准多尺度熵
    mse = None

    # 合成多尺度熵
    cmse = None

    # 精炼的合成多尺度熵
    rcmse = None

    # Correlation Dimension D2
    cd = None

    # Higuchi's Fractal Dimension
    hfd = None

    # Katz's Fractal Dimension
    kfd = None

    # Lempel Ziv Complexity
    lzc = None


class NonLinearFragmentationIndices(object):
    """
    非线性分析中的分段相关指标
    """
    # 序列中拐点数量在整个序列数量的比例
    pip = None
    # 序列中连续加速和减速的段的平均长度的倒数
    ials = None
    # 长度小于3的连续加速和减速片段占全部加速或者减速片段的比例
    pss = None
    # 拐点中长度大于等于4的连续片段占整个的比例
    pas = None


class NonLinearHeartRateAsymmetryIndices(object):
    """
    庞加莱图中非对称心率变异相关指标
    HRA--HRV在时间上呈现出不对称的现象。具体就是心率的加速和减速在长程和短程分析的贡献度不同。
    """
    # 庞加莱图中，特征线上方的点到特征线的距离与除了特征线上的全部点到特征线的距离的比值
    gi = None
    # 庞加莱图中，特征线上方的点的相位角度与除了特征线上的全部点的相位角度的比值
    si = None
    # 庞加莱图中，特征线上方的点的扇形面积之和与除了特征线上的全部点的扇形面积之和的比值
    ai = None
    # 庞加莱图中，特征线下方的点的数量与除了特征线上的全部点的数量的比值
    pi = None
    # the contributions of heart rate decelerations and accelerations
    #             to short-term HRV, respectively (Piskorski,  2011).
    c1_d = None
    c1_a = None
    # short-term variance of contributions of decelerations
    # (prolongations of RR intervals) and accelerations (shortenings of RR intervals),espectively (Piskorski,  2011)
    sd1_d = None
    sd1_a = None
    # the contributions of heart rate decelerations and accelerations to long-term HRV, respectively (Piskorski,  2011).
    c2_d = None
    c2_a = None
    # long-term variance of contributions of decelerations (prolongations of RR intervals)
    # and accelerations (shortenings of RR intervals),respectively (Piskorski,  2011).
    sd2_d = None
    sd2_a = None
    # the total contributions of heart rate decelerations and accelerations to HRV.
    c_d = None
    c_a = None
    # total variance of contributions of decelerations (prolongations of RR intervals)
    # and accelerations (shortenings of RR intervals), respectively (Piskorski,  2011).
    sdnn_d = None
    sdnn_a = None


class NonLinearHeartRateDFAIndices(object):
    """
    非线性分析中去趋势波动指标
    """
    # The monofractal detrended fluctuation analysis of the HR signal corresponding to short-term correlations.
    alpha1 = None

    # The monofractal detrended fluctuation analysis of the HR signal corresponding to long-term correlations .
    alpha2 = None

    #  The multifractal detrended fluctuation analysis of the HR signal corresponding to short-term correlations.
    #  ExpRange is the range of singularity exponents, correspoinding to the width of the singularity spectrum.
    alpha1_ExpRange = None

    # The multifractal detrended fluctuation analysis of the HR signal corresponding to long-term correlations.
    # ExpRange is the range of singularity exponents, correspoinding to the width of the singularity spectrum.
    alpha2_ExpRange = None

    # Multifractal DFA. ExpMean is the mean of singularity exponents.
    alpha1_ExpMean = None

    # Multifractal DFA. ExpMean is the mean of singularity exponents.
    alpha2_ExpMean = None

    # The multifractal detrended fluctuation analysis of the HR signalcorresponding to short-term correlations.
    # DimRange is the range of singularity dimensions, correspoinding to the height of the singularity spectrum.
    alpha1_DimRange = None

    # The multifractal detrended fluctuation analysis of the HR signal corresponding to long-term correlations.
    # DimRange is the range of singularity dimensions, correspoinding to the height of the singularity spectrum.
    alpha2_DimRange = None

    # Multifractal DFA. Dimmean is the mean of singularity dimensions.
    alpha1_DimMean = None

    # Multifractal DFA. Dimmean is the mean of singularity dimensions.
    alpha2_DimMean = None


class NonLinearLorenzIndices(object):
    """
    Lorenz Plot指标
    reference `A Detector for a Chronic Implantable Atrial Tachyarrhythmia Monitor`
    """
    # 区域宽度
    bin_radius = 80
    # 截距
    intercept = bin_radius * math.acos(math.pi / 4)
    # Irregularity Evidence
    irregularity_evidence = 0
    # Density Evidence
    density_evidence = 0
    # Anisotropy Evidence
    anisotropy_evidence = 0
    # PAC Evidence
    pac_evidence = 0
    # Origin Count
    origin_count = 0
    # AF Evidence
    af_evidence = 0
