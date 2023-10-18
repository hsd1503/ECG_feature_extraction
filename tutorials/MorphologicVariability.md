# Morphological Variability（MV）

### 特征说明

* 信号的形态变异性量化；
* 先将指定信号在给定的峰值处切片，然后使用自适应重采样至指定长度，再使用DTW算法计算相邻片段的匹配路径，之后计算片段的变化值，最后使用`Lomb-Scargle`进行谱分析，得到量化结果。


### 使用示例

````python
from features.MorphologicVariability.Extractor import morphological_variability

peaks = []
sampling_rate = 125
ecg = []
print(morphological_variability(ecg, sampling_rate, peaks))

````

### 其他说明

1. 原始论文中选取的是10分钟的片段差异积累量。
2. 本文实现上使用的是相邻差异的作为量化结果，可得到量化结果序列，具体使用可根据实际情况再进行改进。

