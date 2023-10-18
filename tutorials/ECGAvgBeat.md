# Features of ECG Average Beat 

## 特征说明
 针对用户所提供的平均波形信息计算相关的平均波形特征； 特征说明如下：
  
  |  Parameter  | Unit | Description|
  |  ----  | ----  |----|
  | HR        | bpm | 平均心率 |
  | QT        | ms | QT间期 |
  | QTc       | ms | QTc时长 |
  | R/T       | 1 | R波与T波振幅比值 |
  |P Wave|||
  |  exist            | /  | 是否存在 |
  |  form             | /  | P波形态，向上、向下、正负双向、负正双向 |
  |  amplitude        | mV  | P波振幅 |
  |  direction        | /   | 主波方向，向上、向下 |
  |  duration         | ms  | P波时限 |
  |QRS Complex|||
  |  QAmplitude       | mV | Q波振幅 |
  |  RAmplitude       | mV | R波振幅 |
  |  R'Amplitude      | mV | R'波振幅 |
  |  SAmplitude       | mV | S波振幅 |
  |  QDuration        | mV | Q波时限 |
  |  RDuration        | mV | R波时限 |
  |  R'Duration       | mV | R'波时限 |
  |  SDuration        | mV | S波时限 |
  |  QRSPrinciple     | /  | QRS波群主波方向，+向上,-向下 |
  |  QRSAmplitude     | mV | QRS波群振幅 |
  |  QRSDuration      | ms | QRS波群时限 |
  |  Q/R              | 1  | Q波与R波振幅比值 |
  |T Wave|||
  |  exist            | /  | 是否存在 |
  |  form             | /  | T波形态，向上、向下、正负双向、负正双向 |
  |  amplitude        | mV | T波振幅 |
  |  direction        | /  | 主波方向，向上、向下 |
  |  duration         | ms | T波时限 |
  |PR Fragment|||
  |  computable       | /  | 是否可计算 |
  |  interval         | ms | PR间期时长 |
  |  form             | /  | PR段形态，水平、上斜、下斜 |
  |ST Fragment|||
  |  computable       | /  | 是否可计算 |
  |  duration         | ms | ST段时长 |
  |  form             | /  | PR段形态，水平、上斜、下斜 |
  |  amplitude        | mV | ST段上斜、下斜的幅值 |


----------


## 使用样例

````python
from tools.QRSDetector import simple_qrs_detector
from features.ECGAvgBeat.Tools import extract_avg_wave
from features.ECGAvgBeat.FeatureExtractor import extract_features
from tools.SingleBeatBounds import dwt_ecg_delineator

ecg = []
fs = 125

qrs = simple_qrs_detector(ecg, fs)
avg_wave, r = extract_avg_wave(ecg, qrs, fs)

bound_infos = dwt_ecg_delineator(avg_wave, r, fs)

p_start = bound_infos['ECG_P_Onset']
p_end = bound_infos['ECG_P_Offset']
qrs_start = bound_infos['ECG_R_Onset']
qrs_end = bound_infos['ECG_R_Offset']
t_start = bound_infos['ECG_T_Onset']
t_end = bound_infos['ECG_T_Offset']
print(extract_features(avg_wave, p_start, p_end, qrs_start, qrs_end, t_start, t_end, fs))

````
