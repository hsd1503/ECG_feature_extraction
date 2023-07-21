import FrequencyDomain as fd
import TimeDomain as td
import NonLinear as nl
import numpy as np
import time

if __name__ == '__main__':
    peaks = [148, 260, 371, 482, 595, 705, 782, 906, 1016, 1129, 1241,
             1352, 1465, 1576, 1688, 1798, 1911, 2021, 2133, 2244, 2356, 2468,
             2579, 2691, 2803, 2914, 3026, 3137, 3248, 3361, 3472, 3583, 3696,
             ]
    sampling_rate = 125
    rri = np.diff(peaks) / sampling_rate * 1000
    rri = [int(_) for _ in rri]
    print()
    start = time.time()
    print('Frequency Domain Frequencies')
    print(fd.frequencies(peaks, sampling_rate).__dict__)
    print()
    print('Time Domain Basic Statistics')
    print(td.basic_statistics(peaks, sampling_rate).__dict__)
    print()
    print('Time Domain HRV')
    print(td.hrv(peaks, sampling_rate, rri=rri).__dict__)
    print()
    print('Non-Linear Poincare')
    print(nl.poincare(peaks, sampling_rate, rri=rri).__dict__)
    print()
    print('Non-Linear Entropy')
    print(nl.entropy(peaks, sampling_rate, rri=rri).__dict__)
    print()
    print('Non-Linear Fragmentation')
    print(nl.fragmentation(peaks, sampling_rate, rri=rri).__dict__)
    print()
    print('Non-Linear HRA')
    print(nl.hra(peaks, sampling_rate, rri=rri).__dict__)
    print()
    print('Non-Linear DFA')
    print(nl.dfa(peaks, sampling_rate, rri=rri).__dict__)
    print(int((time.time() - start) * 1000))
