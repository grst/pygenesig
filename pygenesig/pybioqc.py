from scipy import stats
import numpy as np


def wmw_u(x, y):
    """
    The Wilcoxon-Mann-Whitney U-statistics.

    Args:
        x (array-like): sample1
        y (array-like): sample2

    Returns:
        float: U

    >>> x = np.arange(0, 10)
    >>> y = np.arange(11, 21)
    >>> u, p = stats.mannwhitneyu(x, y)
    >>> u2 = wmw_u(x, y)
    >>> u == u2
    True
    >>> u
    0.0
    >>> x = np.random.randint(-100, 100, 500)
    >>> y = np.random.randint(-100, 100, 500)
    >>> u, p = stats.mannwhitneyu(x, y)
    >>> u2 = wmw_u(x, y)
    >>> u == u2
    True
    """
    n1 = x.size
    values = np.concatenate((np.array(x), np.array(y)))
    ranks = stats.rankdata(values)
    u = np.sum(ranks[np.arange(0, n1)]) - (n1 * (n1 + 1)) / 2
    return u


def wmw_u_exp(rel_median, n_x, n_y):
    """
    Expected U statistics based on the median rank of a signature.

    Args:
        rel_median:
        n_x:
        n_y:

    Returns:

    """
    median = rel_median * n_y
    return n_x * median - (n_x * (n_x + 1)) / 2


def wmw_r(x, y, u):
    """
    r effect size.

    corrects the BioQC score for size of the signatures.

    Args:
        x:
        y:
        u:

    Returns:

    >>> x = np.arange(0, 10)
    >>> y = np.arange(11, 21)
    >>> u = wmw_u(x, y)
    >>> wmw_r(x, y, u)
    1.0
    >>> u = wmw_u(y, x)
    >>> wmw_r(y, x, u)
    -1.0
    """
    return 1 - (2 * u) / (x.size * y.size)


def wmw_mu(x, y):
    """
    Mean of the normal distribution of U

    Args:
        x:
        y:

    Returns:

    """
    return x.size * y.size / 2


def wmw_sigma(x, y, ranks):
    """
    sigma of the normal distribution of U

    TODO: untested!!

    Args:
        x:
        y:
        ranks:

    Returns:

    """
    n1 = x.size
    n2 = y.size
    n = n1 + n2

    k = 0
    t = []
    r = None
    for rank in ranks:
        if rank != r:
            r = rank
            k += 1
            t.append(1)
        else:
            t[-1] += 1

    sum_corr = 0
    for i in range(k):
        sum_corr += (
            (t[i] ** 3 - t[i]) /
            (n * (n - 1))
        )

    return np.sqrt(
        ((n1 * n2) / 12) * (
            (n + 1) - sum_corr
        )
    )


def mix(x, y, prop_x):
    """
    in-silico mixing of two samples.

    Args:
        x:
        y:
        prop_x:

    Returns:

    """
    return x * prop_x + y * (1 - prop_x)


def wmw_r_corr(x, y, u, u_exp):
    """
    corrected r effect size.

    corrected for signature size and median rank of a signature.

    Args:
        x:
        y:
        u:
        u_exp:

    Returns:

    """
    return 1 - (2 * (u - u_exp)) / (x.size * y.size)