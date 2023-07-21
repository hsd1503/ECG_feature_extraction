from features.HeartRateVariability import FrequencyDomain, TimeDomain, NonLinear


def hrvs(peaks, fs):
    excluded_keys = ['ulf_band', 'vlf_band', 'lf_band', 'hf_band', 'vhf_band', 'psd_method']
    hrv_f = FrequencyDomain.frequencies(peaks, fs).__dict__
    for k in excluded_keys:
        hrv_f[k] = None
    hrv_tb = TimeDomain.basic_statistics(peaks, fs).__dict__
    hrv_t = TimeDomain.hrv(peaks, fs).__dict__
    hrv_nl_p = NonLinear.poincare(peaks, fs).__dict__
    hrv_nl_e = NonLinear.entropy(peaks, fs).__dict__
    return hrv_f, hrv_tb, hrv_t, hrv_nl_p, hrv_nl_e
