"""
计算QT Dispersion
"""


def qtd(qt: list):
    """
    计算QT离散度

    Parameters
    ----------
    qt

    Returns
    -------

    """
    return max(qt) - min(qt)
