"""
非线性分析
"""
import time
from warnings import warn

from features.HeartRateVariability.common import NonLinearPoincareIndices, NonLinearEntropyIndices, \
    NonLinearFragmentationIndices, \
    NonLinearHeartRateAsymmetryIndices, NonLinearHeartRateDFAIndices, hrv_get_rri, NonLinearLorenzIndices
import numpy as np
from neurokit2.signal import signal_zerocrossings
from neurokit2.misc import find_consecutive
from neurokit2.complexity import entropy_approximate, entropy_multiscale, entropy_fuzzy, entropy_sample, \
    entropy_shannon, fractal_correlation, fractal_katz, complexity_lempelziv
from features.HeartRateVariability.fractal_higuchi import fractal_higuchi
from features.HeartRateVariability.fractal_dfa import fractal_dfa


def poincare(peaks: list, sampling_rate=1000, rri: list = None) -> NonLinearPoincareIndices:
    """
    庞加莱相关分析
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `NonLinearPoincareIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]
    pis = NonLinearPoincareIndices()
    rri_n = rrs[:-1]
    rri_plus = rrs[1:]
    x1 = (rri_n - rri_plus) / np.sqrt(2)  # Eq.7
    x2 = (rri_n + rri_plus) / np.sqrt(2)
    pis.sd1 = np.std(x1, ddof=1)
    pis.sd2 = np.std(x2, ddof=1)

    # SD1 / SD2
    pis.sd1_sd2 = pis.sd1 / pis.sd2

    # Area of ellipse described by SD1 and SD2
    pis.s = np.pi * pis.sd1 * pis.sd2

    # CSI / CVI
    T = 4 * pis.sd1
    L = 4 * pis.sd2
    pis.csi = L / T
    pis.cvi = np.log10(L * T)
    pis.csi_modified = L ** 2 / T
    return pis


def entropy(peaks: list, sampling_rate=1000, rri: list = None, **kwargs) -> NonLinearEntropyIndices:
    """
    熵相关分析
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `NonLinearEntropyIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    rrs = np.array(rri)
    rri = rrs[~np.isnan(rrs)]
    # Complexity
    r = 0.2 * np.std(rri, ddof=1)
    ens = NonLinearEntropyIndices()
    try:
        ens.approximate = entropy_approximate(rri, delay=1, dimension=2, tolerance=r)[0]
    except:
        ens.approximate = np.nan

    try:
        ens.sample = entropy_sample(rri, delay=1, dimension=2, tolerance=r)[0]
    except:
        ens.sample = np.nan

    try:
        ens.shannon = entropy_shannon(rri)[0]
    except:
        ens.shannon = np.nan

    try:
        ens.fuzzy = entropy_fuzzy(rri, delay=1, dimension=2, tolerance=r)[0]
    except:
        ens.fuzzy = np.nan

    try:
        ens.mse = entropy_multiscale(rri, dimension=2, tolerance=r, composite=False, refined=False)[0]
    except:
        ens.mse = np.nan

    try:
        ens.cmse = entropy_multiscale(rri, dimension=2, tolerance=r, composite=True, refined=False)[0]
    except:
        ens.cmse = np.nan

    try:
        ens.rcmse = entropy_multiscale(rri, dimension=2, tolerance=r, composite=True, refined=True)[0]
    except:
        ens.rcmse = np.nan

    try:
        ens.cd = fractal_correlation(rri, delay=1, dimension=2, **kwargs)[0]
    except:
        ens.cd = np.nan

    try:
        ens.hfd = fractal_higuchi(rri, **kwargs)
    except:
        ens.hfd = np.nan

    try:
        ens.kfd = fractal_katz(rri)[0]
    except:
        ens.kfd = np.nan

    try:
        ens.lzc = complexity_lempelziv(rri, **kwargs)[0]
    except:
        ens.lzc = np.nan

    return ens


def fragmentation(peaks: list, sampling_rate=1000, rri: list = None) -> NonLinearFragmentationIndices:
    """
    分段相关分析
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `NonLinearFragmentationIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]
    fis = NonLinearFragmentationIndices()
    diff_rri = np.diff(rrs)
    zero_crossings = signal_zerocrossings(diff_rri)

    # Percentage of inflection points (PIP)
    fis.pip = len(zero_crossings) / len(rri)

    # Inverse of the average length of the acceleration/deceleration segments (IALS)
    accelerations = np.where(diff_rri > 0)[0]
    decelerations = np.where(diff_rri < 0)[0]
    consecutive = find_consecutive(accelerations) + find_consecutive(decelerations)
    lengths = [len(i) for i in consecutive]
    fis.ials = 1 / np.average(lengths)

    # Percentage of short segments (PSS) - The complement of the percentage of NN intervals in
    # acceleration and deceleration segments with three or more NN intervals
    fis.pss = np.sum(np.asarray(lengths) < 3) / len(lengths)

    # Percentage of NN intervals in alternation segments (PAS). An alternation segment is a sequence
    # of at least four NN intervals, for which heart rate acceleration changes sign every beat. We note
    # that PAS quantifies the amount of a particular sub-type of fragmentation (alternation). A time
    # series may be highly fragmented and have a small amount of alternation. However, all time series
    # with large amount of alternation are highly fragmented.
    alternations = find_consecutive(zero_crossings)
    lengths = [len(i) for i in alternations]
    fis.pas = np.sum(np.asarray(lengths) >= 4) / len(lengths)
    return fis


def hra(peaks: list, sampling_rate=1000, rri: list = None) -> NonLinearHeartRateAsymmetryIndices:
    """
    心率不对称相关分析
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `NonLinearHeartRateAsymmetryIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]
    N = len(rrs) - 1
    # rri_n, x-axis
    x = rrs[:-1]
    # rri_plus, y-axis
    y = rrs[1:]

    diff = y - x
    # set of points above IL where y > x
    decelerate_indices = np.where(diff > 0)[0]
    # set of points below IL where y < x
    accelerate_indices = np.where(diff < 0)[0]
    nochange_indices = np.where(diff == 0)[0]

    # Distances to centroid line l2
    centroid_x = np.mean(x)
    centroid_y = np.mean(y)
    dist_l2_all = abs((x - centroid_x) + (y - centroid_y)) / np.sqrt(2)

    # Distances to LI
    dist_all = abs(y - x) / np.sqrt(2)

    # Calculate the angles
    # phase angle LI - phase angle of i-th point
    theta_all = abs(np.arctan(1) - np.arctan(y / x))
    # Calculate the radius
    r = np.sqrt(x ** 2 + y ** 2)
    # Sector areas
    S_all = 1 / 2 * theta_all * r ** 2

    # Guzik's Index (GI)
    den_GI = np.sum(dist_all)
    num_GI = np.sum(dist_all[decelerate_indices])

    hras = NonLinearHeartRateAsymmetryIndices()
    hras.gi = (num_GI / den_GI) * 100

    # Slope Index (SI)
    den_SI = np.sum(theta_all)
    num_SI = np.sum(theta_all[decelerate_indices])
    hras.si = (num_SI / den_SI) * 100

    # Area Index (AI)
    den_AI = np.sum(S_all)
    num_AI = np.sum(S_all[decelerate_indices])
    hras.ai = (num_AI / den_AI) * 100

    # Porta's Index (PI)
    m = N - len(nochange_indices)  # all points except those on LI
    b = len(accelerate_indices)  # number of points below LI
    if m == 0:
        hras.pi = np.nan
    else:
        hras.pi = (b / m) * 100

        # Short-term asymmetry (SD1)
    sd1d = np.sqrt(np.sum(dist_all[decelerate_indices] ** 2) / (N - 1))
    sd1a = np.sqrt(np.sum(dist_all[accelerate_indices] ** 2) / (N - 1))

    sd1I = np.sqrt(sd1d ** 2 + sd1a ** 2)
    hras.c1_d = (sd1d / sd1I) ** 2
    hras.c1_a = (sd1a / sd1I) ** 2
    hras.sd1_d = sd1d  # SD1 deceleration
    hras.sd1_a = sd1a  # SD1 acceleration
    # out["SD1I"] = sd1I  # SD1 based on LI, whereas SD1 is based on centroid line l1

    # Long-term asymmetry (SD2)
    longterm_dec = np.sum(dist_l2_all[decelerate_indices] ** 2) / (N - 1)
    longterm_acc = np.sum(dist_l2_all[accelerate_indices] ** 2) / (N - 1)
    longterm_nodiff = np.sum(dist_l2_all[nochange_indices] ** 2) / (N - 1)

    sd2d = np.sqrt(longterm_dec + 0.5 * longterm_nodiff)
    sd2a = np.sqrt(longterm_acc + 0.5 * longterm_nodiff)

    sd2I = np.sqrt(sd2d ** 2 + sd2a ** 2)
    hras.c2_d = (sd2d / sd2I) ** 2
    hras.c2_a = (sd2a / sd2I) ** 2
    hras.sd2_d = sd2d  # SD2 deceleration
    hras.sd2_a = sd2a  # SD2 acceleration
    # out["SD2I"] = sd2I  # identical with SD2

    # Total asymmerty (SDNN)
    sdnnd = np.sqrt(0.5 * (sd1d ** 2 + sd2d ** 2))  # SDNN deceleration
    sdnna = np.sqrt(0.5 * (sd1a ** 2 + sd2a ** 2))  # SDNN acceleration
    sdnn = np.sqrt(sdnnd ** 2 + sdnna ** 2)  # should be similar to sdnn in hrv_time
    hras.c_d = (sdnnd / sdnn) ** 2
    hras.c_a = (sdnna / sdnn) ** 2
    hras.sdnn_d = sdnnd
    hras.sdnn_a = sdnna
    return hras


def dfa(peaks: list, sampling_rate=1000, rri: list = None, n_windows="default",
        **kwargs) -> NonLinearHeartRateDFAIndices:
    """
    去趋势波动相关分析
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param n_windows: 窗口宽度
    :param rri: rr间隙,单位ms
    :return: see `NonLinearHeartRateDFAIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]

    if 'dfa_windows' in kwargs:
        dfa_windows = kwargs['dfa_windows']
    else:
        dfa_windows = [(4, 11), (12, None)]

    # Determine max beats
    if dfa_windows[1][1] is None:
        max_beats = len(peaks) / 10
    else:
        max_beats = dfa_windows[1][1]

    # No. of windows to compute for short and long term
    if n_windows == "default":
        n_windows_short = int(dfa_windows[0][1] - dfa_windows[0][0] + 1)
        n_windows_long = int(max_beats - dfa_windows[1][0] + 1)
    elif isinstance(n_windows, list):
        n_windows_short = n_windows[0]
        n_windows_long = n_windows[1]
    dfas = NonLinearHeartRateDFAIndices()
    # Compute DFA alpha1
    short_window = np.linspace(dfa_windows[0][0], dfa_windows[0][1], n_windows_short).astype(int)
    # For monofractal
    d = fractal_dfa(rrs, multifractal=False, windows=short_window, **kwargs)
    if isinstance(d, float) and np.isnan(d):
        dfas.alpha1 = np.nan
    else:
        dfas.alpha1 = d['slopes'][0]
    # For multifractal
    mdfa_alpha1 = fractal_dfa(rrs, multifractal=True, q=np.arange(-5, 6), windows=short_window, **kwargs)
    if isinstance(mdfa_alpha1, float) and np.isnan(mdfa_alpha1):
        dfas.alpha1_ExpRange = np.nan
        dfas.alpha1_ExpMean = np.nan
        dfas.alpha1_DimRange = np.nan
        dfas.alpha1_DimMean = np.nan
    else:
        dfas.alpha1_ExpRange = mdfa_alpha1['ExpRange']
        dfas.alpha1_ExpMean = mdfa_alpha1['ExpMean']
        dfas.alpha1_DimRange = mdfa_alpha1['DimRange']
        dfas.alpha1_DimMean = mdfa_alpha1['DimMean']

    # Compute DFA alpha2
    # sanatize max_beats
    if max_beats < dfa_windows[1][0] + 1:
        warn(
            "DFA_alpha2 related indices will not be calculated. "
            "The maximum duration of the windows provided for the long-term correlation is smaller "
            "than the minimum duration of windows. Refer to the `windows` argument in `nk.fractal_dfa()` "
            "for more information."
        )
    else:
        long_window = np.linspace(dfa_windows[1][0], int(max_beats), n_windows_long).astype(int)
        # For monofractal
        dfas.alpha2 = fractal_dfa(rrs, multifractal=False, windows=long_window, **kwargs)['slopes'][0]
        # For multifractal
        mdfa_alpha2 = fractal_dfa(rrs, multifractal=True, q=np.arange(-5, 6), windows=long_window, **kwargs)
        dfas.alpha2_ExpRange = mdfa_alpha2['ExpRange']
        dfas.DFA_alpha2_ExpMean = mdfa_alpha2['ExpMean']
        dfas.DFA_alpha2_DimRange = mdfa_alpha2['DimRange']
        dfas.alpha2_DimMean = mdfa_alpha2['DimMean']
    return dfas


def lorenz(peaks: list, sampling_rate=1000, rri: list = None) -> NonLinearLorenzIndices:
    """
    Lorenz相关分析
    :param peaks: R波下标位置
    :param sampling_rate: 采样率
    :param rri: rr间隙,单位ms
    :return: see `NonLinearLorenzIndices`
    """
    if rri is None:
        rri, sampling_rate = hrv_get_rri(peaks, sampling_rate)
    rrs = np.array(rri)
    rrs = rrs[~np.isnan(rrs)]
    deltas = np.diff(rrs)
    lorens = NonLinearLorenzIndices()
    point_counts = np.zeros(13)
    for idx in range(1, len(deltas)):
        x = deltas[idx]
        y = deltas[idx - 1]
        if -lorens.bin_radius < x < lorens.bin_radius and -lorens.bin_radius < y < lorens.bin_radius:
            point_counts[0] += 1
        elif x * y == 0:
            if x != 0:
                if x > 0:
                    point_counts[7] += 0
                else:
                    point_counts[5] += 0
            else:
                if y > 0:
                    point_counts[8] += 0
                else:
                    point_counts[6] += 0
        else:
            if x * y > 0:
                if x > 0:
                    if x + lorens.intercept < y:
                        if x > lorens.bin_radius:
                            point_counts[4] += 1
                        else:
                            point_counts[8] += 1
                    elif x - lorens.intercept > y:
                        if y > lorens.bin_radius:
                            point_counts[4] += 1
                        else:
                            point_counts[7] += 1
                    else:
                        point_counts[12] += 1
                else:
                    if x + lorens.intercept < y:
                        if y > -lorens.bin_radius:
                            point_counts[5] += 1
                        else:
                            point_counts[2] += 1
                    elif x - lorens.intercept > y:
                        if x > -lorens.bin_radius:
                            point_counts[6] += 1
                        else:
                            point_counts[2] += 1
                    else:
                        point_counts[10] += 1
            else:
                if x > 0:
                    if -x + lorens.intercept < y:
                        if y > -lorens.bin_radius:
                            point_counts[7] += 1
                        else:
                            point_counts[3] += 1
                    elif -x - lorens.intercept > y:
                        if x > lorens.bin_radius:
                            point_counts[3] += 1
                        else:
                            point_counts[6] += 1
                    else:
                        point_counts[11] += 1

                else:
                    if -x + lorens.intercept < y:
                        if x < -lorens.bin_radius:
                            point_counts[1] += 1
                        else:
                            point_counts[8] += 1
                    elif -x - lorens.intercept > y:
                        if y > lorens.bin_radius:
                            point_counts[1] += 1
                        else:
                            point_counts[5] += 1
                    else:
                        point_counts[9] += 1

    lorens.irregularity_evidence = len(np.where(point_counts[1:] > 0)[0])
    lorens.density_evidence = sum([_ - (0 if _ < 1 else 1) for _ in point_counts[5:]])
    lorens.anisotropy_evidence = abs(sum(point_counts[[9, 11]]) - sum(point_counts[[10, 12]])) \
                                 + abs(sum(point_counts[[6, 7]]) - sum(point_counts[[5, 8]]))
    lorens.pac_evidence = sum([_ - (0 if _ < 1 else 1) for _ in point_counts[1:5]]) \
                          + sum([_ - (0 if _ < 1 else 1) for _ in point_counts[[5, 6, 10]]]) \
                          - sum([_ - (0 if _ < 1 else 1) for _ in point_counts[[7, 8, 12]]])
    lorens.origin_count = len(deltas) - 1
    lorens.af_evidence = lorens.irregularity_evidence - lorens.origin_count - 2 * lorens.pac_evidence
    return lorens


if __name__ == '__main__':
    import tftb.processing as tfp
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from neurokit2 import rsp

    import neurokit2 as nk
    import mne

    # Download data
    data = nk.data("bio_resting_5min_100hz")
    # Process the data
    df, info = nk.rsp_process(data["RSP"], sampling_rate=100)

    # Single dataframe is passed
    nk.rsp_intervalrelated(df, sampling_rate=100) #doctest: +SKIP

    epochs = nk.epochs_create(df, events=[0, 15000], sampling_rate=100, epochs_end=150)
    nk.rsp_intervalrelated(epochs) #doctest: +SKIP
    print()

    # f, axx = plt.subplots(1, 1)
    #
    # x = np.array(
    #     [94, 153, 209, 277, 345, 434, 488, 583, 652, 712, 789, 846, 951, 1028, 1098, 1206, 1275, 1339, 1412, 1493, 1700,
    #      1781, 1862, 1945, 2020, 2083, 2186, 2246, 2317, 2377, 2479, 2579, 2634, 2698, 2785, 2898, 2958, 3025, 3097,
    #      3206, 3256, 3345, 3408, 3461, 3555, 3640, 3702, 3759, 3822, 3917, 3980, 4050, 4150, 4211, 4265, 4378, 4459,
    #      4562, 4638, 4691, 4777, 4859, 4942, 5016, 5113, 5204, 5266, 5323, 5374, 5432, 5491, 5550, 5645, 5739, 5836,
    #      5935, 6023, 6114, 6178, 6248, 6317, 6419, 6493, 6559, 6613, 6717, 6778, 6825, 6865])
    # x, f = hrv_get_rri(x.tolist(), sampling_rate=125, interpolate=True)
    # # x *= (1 / (np.max(x) - np.min(x)))
    # # x = pd.read_csv('c:/Users/Desktop/demo.txt', header=None).values[0]
    # start = time.time()
    # transformed = tfp.smoothed_pseudo_wigner_ville(np.array(x))
    # print('耗时{}ms'.format(int(time.time() * 1000 - start * 1000)))
    # sns.set(font_scale=1.5)
    # sns.set_context({"figure.figsize": (8, 8)})
    # sns.heatmap(data=pd.DataFrame(np.array(transformed)), square=True)
    # plt.show()
    # print()
