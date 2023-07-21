# Features of QRS Axis

## 特征说明

通过用户给定的多导联心电参数，计算QRS电轴，[参考链接](https://en.my-ekg.com/calculation-ekg/heart-axis-calculator.php) 。

----------

## 使用样例

````python
from features.QRSAxis.FeatureExtractor import *

q_amplitudes = [1, 2]
r_amplitudes = [4, 6]
s_amplitudes = [None, 0]
q_r_s_amplitude_i = calculate_q_r_s_amplitude(q_amplitudes[0], r_amplitudes[0], s_amplitudes[0])
q_r_s_amplitude_iii = calculate_q_r_s_amplitude(q_amplitudes[1], r_amplitudes[1], s_amplitudes[1])
print(calculate_axis_using_i_iii_lead(q_r_s_amplitude_i, q_r_s_amplitude_iii))

````
