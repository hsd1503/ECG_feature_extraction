# Heart Rate Variability（HRV）

### 特征说明

* 健康生物的系统模式是复杂且可变的，可以用数理上的混乱度来描述。心率变异性（HRV）包括连续心跳之间的时间间隔的变化，称为心跳间期（IBI）。
* 一个健康的心脏不是节拍器。健康心脏的振荡是复杂且不断变化的，它使心血管系统能够迅速适应来自生理和心理的突然刺激。
* 对监测期间观察到的心率变异的数量进行量化，监测范围从2分钟到24小时不等。频域值计算特定带宽内信号能量的绝对或相对数量。非线性测量对一系列IBI的不可预测性和复杂性进行量化。
* 心率变异性测量强调测量环境的重要性，包括记录时长、受试者的年龄和性别等，对心率变异性基础值的影响。
* 可以用时域、频域和非线性测量来描述24小时、短期（ST，约5分钟）或短暂以及超短期（UST，<5分钟）心率变异。
  由于较长记录时间能更好反应波动较慢的过程（如昼夜节律）和心血管系统对更丰富的环境刺激和工作负荷的反应，所以短期和超短期的数值不能与24小时的数值互换。

--------

* Time-Domain Measures


  |  Parameter  | Unit | Description|
  |  ----  | ----  |----|
  |Basic Statistics|||
  | max_rr        | ms | 最大RR间隙 |
  | min_rr        | ms | 最小RR间隙 |
  | mean_rr       | ms | 平均RR间隙 |
  | median_rr     | ms | RR间隙中位数 |
  | mad_rr        | ms | 中值绝对偏差 |
  |  mcv_rr       | ms | 归一化后的中值绝对偏差 |
  |  variance_rr  | ms | RR间隙方差 |
  |  std_rr       | ms | RR间隙母体标准差（自由度为n） |
  |  amplitude_rr | ms | RR间隙波动范围 |
  | percentile_5  | ms | RR间隙5%分位数 |
  | percentile_25 | ms | RR间隙25%分位数 |
  | percentile_75 | ms | RR间隙75%分位数 |
  | percentile_95 | ms | RR间隙95%分位数 |
  | percentile_range_75_25| ms | RR间隙75%分位数与25%分位数差值|
  | percentile_range_95_5 | ms | RR间隙95%分位数与5%分位数差值|
  | skew_rr       | ms | RR间隙序列偏度 |
  | kurtosis_rr   | ms | RR间隙序列峰度 |
  | count_rr      | /  | RR间隙数量 |
  | Standard Indices|||
  | sdnn          | ms | RR间隙样本标准差（自由度为n-1）|
  | sdann_1       | ms | 连续的1分钟划分一段，每段的RR间隙均值的样本标准差 |
  | sdann_2       | ms | 连续的2分钟划分一段，每段的RR间隙均值的样本标准差 |
  | sdann_5       | ms | 连续的5分钟划分一段，每段的RR间隙均值的样本标准差 |
  | sdnni_1       | ms | 连续的1分钟划分一段，每段的RR间隙sdnn的均值 |
  | sdnni_2       | ms | 连续的2分钟划分一段，每段的RR间隙sdnn的均值 |
  | sdnni_5       | ms | 连续的5分钟划分一段，每段的RR间隙sdnn的均值 |
  | rmssd         | ms | 连续RR间隙差值的均方根 |
  | sdsd          | ms | 连续RR间隙差值的样本标准差 |
  | cvnn          | 1  | sdnn与mean_rr比值 |
  | cvsd          | 1  | rmssd与mean_rr比值 |
  | pNN50         | %  | 计算连续RR间隙差值中，大于等于50ms的数量百分比 |
  | pNN20         | %  | 计算连续RR间隙差值中，大于等于20ms的数量百分比 |
  | tri_index     | 1  | RR间隙数量与RR间隙直方图中最大值的比值  |
  | tinn          | ms | 拟合RR间隙直方图中最佳三角形对应的底边差值 |
  | tinn_m        | ms | tinn的底边端点m的值 |
  | tinn_n        | ms | tinn的底边端点n的值 |

-------

* Frequency-Domain Measures


  |  Parameter  | Unit | Description|
  |  ----       | ---- |----        |
  | ulf_power   | ms²  | 超低频率（0-0.0033Hz）功率密度积分|
  | vlf_power   | ms²  | 较低频率（0.0033-0.04Hz）功率密度积分|
  | lf_power    | ms²  | 低频率（0.04-0.15Hz）功率密度积分 |
  | hf_power    | ms²  | 高频率（0.15-0.4Hz）功率密度积分|
  | vhf_power   | ms²  | 较高频率（0.4-0.5Hz）功率密度积分 |
  | total_power | ms²  | 0-0.5Hz频带内的功率密度积分 |
  | lf_n        | 1    | lf_power与total_power比值 |
  | hf_n        | 1    | hf_power与total_power比值 |
  | lf_hf       | 1    | lf_power与hf_power比值  |

-----------

* Non-Linear Measures


|  Parameter  | Unit | Description|
|  ----       | ---- |----        |
| Poincare Indices|  |            |
| s           | ms²  | 庞加莱图中椭圆面积 |
| sd1         | ms   | 庞加莱图中的椭圆上垂直于特征线最大的线段长度|
| sd2         | ms   | 庞加莱图中的椭圆上平行于特征线最大的线段长度|
| sd1_sd2     | 1    | sd1与sd2比值 |
| csi         | 1    | 通过庞加莱图计算出的交感神经指数 |
| cvi         | ms²  | 通过庞加莱图计算出的迷走神经指数 |
| csi_modified| ms   | 修正后的交感神经指数 |
| Entropy Indices|||
| approximate | /    | RR间隙近似熵 |
| sample      | /    | RR间隙样本熵 |
| shannon     | /    | RR间隙香侬熵 | 
| fuzzy       | /    | RR间隙模糊熵 |
| mse         | /    | RR间隙标准多尺度熵 |
| cmse        | /    | RR间隙合成多尺度熵 |
| rcmse       | /    | RR间隙精炼的合成多尺度熵 |
| cd          | /    | Correlation Dimension D2 |
| hfd         | /    | Higuchi's Fractal Dimension |
| kfd         | /    | Katz's Fractal Dimension |
| lzc         | /    | Lempel Ziv Complexity|
| Fragmentation Indices|||
| pip         | 1    | 序列中拐点数量在整个序列数量的比例 |
| ials        | 1    | 序列中连续加速和减速的段的平均长度的倒数 |
| pss         | 1    | 长度小于3的连续加速和减速片段占全部加速或者减速片段的比例 |
| pas         | 1    | 拐点中长度大于等于4的连续片段占整个的比例 |
| Heart Rate Asymmetry Indices|||
| gi          | /    | 庞加莱图中，特征线上方的点到特征线的距离与除了特征线上的全部点到特征线的距离的比值 |
| si          | /    | 庞加莱图中，特征线上方的点的相位角度与除了特征线上的全部点的相位角度的比值 |
| ai          | /    |  庞加莱图中，特征线上方的点的扇形面积之和与除了特征线上的全部点的扇形面积之和的比值|
| pi         | /    |  庞加莱图中，特征线下方的点的数量与除了特征线上的全部点的数量的比值|
| c1_d       | /   |  减速在短期HRV中的贡献度 |
| c1_a       | /   |  加速在短期HRV中的贡献度 |
| sd1_d      | /   |  减速在短期HRV方差中的贡献度 |
| sd1_a      | /   |  加速在短期HRV方差中的贡献度 |
| c2_d       | /   |  减速在长期HRV中的贡献度 |
| c2_a       | /   |  加速在长期HRV中的贡献度 |
| sd2_d      | /   |  减速在长期HRV方差中的贡献度 |
| sd2_a      | /   |  加速在长期HRV方差中的贡献度 |
| c_d        | /   |  减速在HRV中的总贡献度 |
| c_a        | /   |  加速在HRV中的总贡献度 |
| sdnn_d     | /   |  减速在总体方差中的贡献度 |
| sdnn_a     | /   |  加速在总体方差中的贡献度 |
| Detrended Fluctuation Analysis Indices|||
| alpha1     | /   |  单分形去趋势波动分析在短期HRV上的相关性|
| alpha2     | /   |  单分形去趋势波动分析在长期HRV上的相关性|
| alpha1_ExpRange |/| alpha1对应奇异谱的宽度|
| alpha2_ExpRange | / | alpha2对应奇异谱的宽度|
| alpha1_ExpMean | / |  alpha1的均值 |
| alpha2_ExpMean | / |  alpha2的均值 |
| alpha1_DimRange| / | alpha1对应奇异谱的高度|
| alpha2_DimRange| / | alpha2对应奇异谱的高度|
| alpha1_DimMean | / | alpha1奇异点维度的均值|
| alpha2_DimMean | / | alpha2奇异点维度的均值|
|Lorenz Indices|||
|irregularity_evidence|1|由Lorenz图获取的不规则指标|
|density_evidence|1|由Lorenz图获取的密度指标|
|anisotropy_evidence|1|由Lorenz图获取的各向异性指标|
|pac_evidence|1|由Lorenz图获取的早搏指标|
|origin_count|1|Lorenz图中的总统计数量|
|af_evidence|1|由Lorenz图获取的房颤指标|

### 使用示例

````python
import features.HeartRateVariability.FrequencyDomain as fd
import features.HeartRateVariability.TimeDomain as td
import features.HeartRateVariability.NonLinear as nl

peaks = [148, 260, 371, 482, 595, 705, 782, 906, 1016, 1129, 1241,
         1352, 1465, 1576, 1688, 1798, 1911, 2021, 2133, 2244, 2356, 2468,
         2579, 2691, 2803, 2914, 3026, 3137, 3248, 3361, 3472, 3583, 3696]
sampling_rate = 125
print()
print('Frequency Domain Frequencies')
print(fd.frequencies(peaks, sampling_rate).__dict__)
print()
print('Time Domain Basic Statistics')
print(td.basic_statistics(peaks, sampling_rate).__dict__)
print()
print('Time Domain HRV')
print(td.hrv(peaks, sampling_rate).__dict__)
print()
print('Non-Linear Poincare')
print(nl.poincare(peaks, sampling_rate).__dict__)
print()
print('Non-Linear Entropy')
print(nl.entropy(peaks, sampling_rate).__dict__)

# Frequency Domain Frequencies
# {'ulf_band': (0, 0.0033), 'vlf_band': (0.0033, 0.04), 'lf_band': (0.04, 0.15), 'hf_band': (0.15, 0.4), 'vhf_band': (0.4, 0.5), 'psd_method': 'welch', 'total_power': 0.23491727179776173, 'ulf_power': 0, 'vlf_power': 0, 'lf_power': 0, 'hf_power': 0.1628652155447752, 'vhf_power': 0.07205205625298652, 'lf_n': 0.0, 'hf_n': 0.6932875318124095, 'lf_hf': 0.0}

# Time Domain Basic Statistics
# {'max_rr': 992.0, 'min_rr': 616.0, 'mean_rr': 887.0, 'median_rr': 892.0, 'mad_rr': 5.9304, 'mcv_rr': 0.006648430493273542, 'count_rr': 32, 'variance_rr': 2723.0, 'std_rr': 52.182372502598994, 'amplitude_rr': 376.0, 'percentile_5': 880.0, 'percentile_25': 888.0, 'percentile_75': 896.0, 'percentile_95': 904.0, 'percentile_range_75_25': 8.0, 'percentile_range_95_5': 24.0, 'skew_rr': -4.114880005464073, 'kurtosis_rr': 20.24638607830441}

# Time Domain HRV
# {'sdnn': 53.017343480010574, 'sdann_1': nan, 'sdnni_1': nan, 'sdann_2': nan, 'sdnni_2': nan, 'sdann_5': nan, 'sdnni_5': nan, 'rmssd': 85.88664697746826, 'sdsd': 87.30596304259922, 'cvnn': 0.05977152590756547, 'cvsd': 0.09682823785509387, 'pNN20': 22.58064516129032, 'pNN50': 9.67741935483871, 'tri_index': 2.909090909090909, 'tinn_m': 906.25, 'tinn_n': 617.1875, 'tinn': 289.0625}

# Non-Linear Poincare
# {'sd1': 61.734638505443996, 'sd2': 44.55524539796768, 'sd1_sd2': 1.3855750979268522, 's': 8641.27093613642, 'csi': 0.7217219777522246, 'cvi': 4.643547732065502, 'csi_modified': 128.62599931142776}

# Non-Linear Entropy
# {'approximate': 0.3020827959592527, 'sample': 0.29725152346793177, 'shannon': 2.184598505552562, 'fuzzy': 0.4518949863871905, 'mse': nan, 'cmse': nan, 'rcmse': nan, 'cd': 0.17076614634791326, 'kfd': 1.6360218596591611, 'lzc': 1.25}

````

### 其他说明

1. 理论上由于HRV在时间上的对称性，即心率的加速和减速在长程和短程分析的贡献度不同，所以长程和短程分析指标一般不能混用。
2. 部分分析指标使用的是RR间隙，有的使用的是NN间隙，理论上二者亦不能混用。
