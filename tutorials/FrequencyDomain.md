# Frequency-Domain Features

### 特征说明

* 用于提取ECG数据的功率谱信息；
* 输入包括；一维ECG数据及其采样率；
* 输出为（2，3751）第一行是频谱图对应的频率值，第二行是每个频率值所对应的功率值。

### 使用示例

````python
from features.ECGFrequencyDomain.Extractor import psd

ecg = [...]
fs = 125
psd(ecg, fs)
````

