import numpy as np


def smooth_avg1(in_put, radius):
    """
    均值滤波2
    中心值减去其前后半径内的均值作为中心值
    :param in_put: 平滑前的数据,一维数组
    :param radius: 平滑半径
    :return:
    """
    output = []
    # 总循环次数
    size = len(in_put)
    # 没有滤波必要
    if size <= 2 * radius + 1:
        return np.array(in_put)

    # 累加数
    sum_ = 0

    # 分段计算，减少循环判断次数

    # 计算前radius个数据和
    for i in range(radius):
        sum_ += in_put[i]

    # 前radius个数据仅需要加最末尾数据
    for i in range(radius + 1):
        sum_ += in_put[i + radius]
        output.append(in_put[i] - sum_ / (i + radius + 1))

    # radius至size-radius需要
    # 减去当前数的前radius位置的数
    # 再加上当前数的后radius位置的数
    for i in range(radius + 1, size - radius):
        sum_ -= in_put[i - radius - 1]
        sum_ += in_put[i + radius]
        output.append(in_put[i] - sum_ / (2 * radius + 1))

    # 序列末尾的pointsCount个数据仅需要减去其最前面的数据即可
    for i in range(size - radius, size):
        sum_ -= in_put[i - radius - 1]
        output.append(in_put[i] - sum_ / (size - (i - radius)))
    return np.array(output)
