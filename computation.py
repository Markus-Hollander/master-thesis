from fastcache import clru_cache
import numpy as np
from scipy.stats import hypergeom, chi2_contingency, pearsonr, norm
from typing import Any
from util import has_zero, remove_zero


def list_stats(val: list[float], r=2) -> tuple[int, int, float, float, float, float, float]:
    """
    :param val: numpy array with values
    :param r: integer specifying the rounding precision
    :return: n, number non-zero, % non-zero, minimum, median, mean and maximum of the given values
    """
    if len(val) > 0:
        non_zero = np.count_nonzero(val)                        # type: int
        stats = ((non_zero / len(val)) * 100, np.amin(val),
                 np.median(val), np.mean(val), np.amax(val))    # type: tuple[float, float, float, float, float]
        stats = tuple([round(x, r) for x in stats])             # type: tuple[float, float, float, float, float]

        return (len(val), non_zero) + stats

    return 0, 0, 0.0, np.NaN, np.NaN, np.NaN, np.NaN


def benjamini_help(sorted_p: list[float], q: float) -> list[bool]:
    """
    One stage Benjamini-Hochberg procedure to identify all p-values that are statistically significant when controlling
    the false-discovery rate (FDR) in multiple hypothesis testing.
    :param sorted_p: list of p-values in ascending order
    :param q: desired, FDR controlled significance threshold
    :return: significance assessment for each input p-value
    """
    # number of p-values
    n = len(sorted_p)
    # index of the highest p-value that meets the criterion
    k = -1

    # identify the highest p-value that meets the criterion
    for i, p in reversed(list(enumerate(sorted_p))):
        # enumerate index is 0 based
        if p <= ((i + 1) * q) / n:
            k = i
            break

    # label all smaller p-values as statistically significant
    return [True if i <= k else False for i in range(n)]


def benjamini(sorted_p: list[float], q: float, krieger=False) -> list[bool]:
    """
    Performs multiple hypothesis testing correction that controls the false-discovery rate (FDR).
    :param sorted_p: sorted_p: list of p-values in ascending order
    :param q: desired, FDR controlled significance threshold
    :param krieger: True if the two-stage version by Benjamini-Krieger-Yekutieli should be performed,
                    False if one the one stage Benjamini-Hochberg procedure should be performed
    :return: significance assessment for each input p-value
    """
    # number of p-values
    n = len(sorted_p)
    # adjusted the significance threshold for Benjamini-Krieger-Yekutieli
    q_prime = q / (1 + q) if krieger else q
    # one-stage Benjamini-Hochberg to identify statistically significant p-values
    significant = benjamini_help(sorted_p, q_prime)

    # second stage for Benjamini-Krieger-Yekutieli
    if krieger:
        # number of significant p-values in the first stage
        n_sig = significant.count(True)
        # if none were significant, the procedure stops here
        if n_sig == 0:
            return [False] * n
        # if all were significant, the procedure stops here
        if n_sig == n:
            return [True] * n

        # estimated number of true null-hypotheses
        n_null = n - n_sig
        # adjusted significance level
        q_star = (q_prime * n) / n_null
        # second iteration of the one-stage Benjamini-Hochberg procedure
        significant = benjamini_help(sorted_p, q_star)

    return significant


def benjamini_dict(p_values: dict[Any, float], q: float, krieger=False) -> dict[Any, bool]:
    """
    Performs multiple hypothesis testing correction that controls the false-discovery rate (FDR) for p-values
    dictionaries.
    :param p_values: dictionary with p-values as values and some sort of identifiers as keys
    :param q: desired, FDR controlled significance threshold
    :param krieger: True if the two-stage version by Benjamini-Krieger-Yekutieli should be performed,
                    False if one the one stage Benjamini-Hochberg procedure should be performed
    :return: significance assessment for each input p-value
    """
    # sorted key, p-value pairs by the p-value in ascending order and skip pairs where the p-value is not a number
    sorted_p = sorted([(k, p) for k, p in p_values.items() if not is_nan(p)], key=lambda x: x[1])

    # if there were no valid p-values, return False for all input keys
    if not sorted_p:
        return {k: False for k in p_values.keys()}

    # split the keys and p-values in the same order
    keys, sorted_p = zip(*sorted_p)

    # initialise the return dictionary
    ret = {k: False for k in p_values.keys()}

    # perform multiple hypothesis testing correction and match the results to the keys
    for k, p in zip(keys, benjamini(sorted_p, q, krieger)):
        ret[k] = p

    return ret


@clru_cache(maxsize=None)
def hyper_geometric(population: int, population_successes: int, draws: int, draw_successes: int) -> float:
    """
    hypergeom.sf(k - 1, M, n, N)
    :param population: M, population size
    :param population_successes: n, successes in the population
    :param draws: N, number of draws
    :param draw_successes: k, successes in the draw
    :return:
    """
    return hypergeom.sf(k=draw_successes - 1, M=population, n=population_successes, N=draws)


def chi_square(min_n: int, value_pairs: list[tuple[int, int]]) -> tuple[float, float, float, str]:
    """
    Pearson's chi-squared test of independence for two variables with at most 3 categories each. Computes the
    probability of observing the input value pairs under the assumption that they come from independent distributions.
    :param min_n: minimum number of input value pairs
    :param value_pairs: list of categorical value pairs
    :return: p-value, Cramer's v, corrected Cramer's v, comment
    """
    # sample size
    n = len(value_pairs)
    # initialise row and column numbers
    r = 3
    c = 3
    # initialise the outcome comment
    comment = 'success'

    # make sure the sample size is sufficient
    if n < min_n:
        return np.nan, np.nan, np.nan, 'less than {0} values'.format(min_n)

    # initialise the contingency table
    table = [[0 for _ in range(3)] for _ in range(3)]

    # fill the contingency table by counting pairs
    for i, j in value_pairs:
        table[j + 1][i + 1] += 1

    # remove 0 columns and/or 0 rows
    if has_zero(table):
        table = remove_zero(table)

        if len(table) < 1:
            return np.nan, np.nan, np.nan, 'no rows after removing empty rows/columns'

        if len(table[0]) < 1:
            return np.nan, np.nan, np.nan, 'no columns after removing empty rows/columns'

        comment += ' - removed empty rows/columns'

    # based on the contingency table, compute the chi-squared statistic, p-value, degrees of freedom and expected counts
    chi2, p, dof, expected = chi2_contingency(table)

    # compute Cramer's v
    phi_2 = chi2 / n
    v = np.sqrt(phi_2 / min(r - 1, c - 1))

    estimated_phi_2 = max(0.0, phi_2 - (((c - 1) * (r - 1)) / (n - 1)))
    estimated_r = r - (pow(r - 1, 2) / (n - 1))
    estimated_c = c - (pow(c - 1, 2) / (n - 1))

    # additionally compute the corrected Cramer's v
    corrected_v = np.sqrt(estimated_phi_2 / min(estimated_r - 1, estimated_c - 1))

    return p, v, corrected_v, comment


def correlation(min_n: int, value_pairs: list[tuple[int, int]]) -> tuple[float, str]:
    """
    Pearson's correlation coefficient for categorical value pairs.
    :param min_n: minimum number of input value pairs
    :param value_pairs: list of categorical value pairs
    :return: Pearson's correlation coefficient
    """
    # ensure that the sample size is sufficient
    if len(value_pairs) < min_n:
        return np.nan, 'less than {0} values ({1})'.format(min_n, len(value_pairs))

    # obtain a matching list of each variable
    x, y = zip(*value_pairs)

    # compute the unique values for each variable
    x_unique = set(x)
    y_unique = set(y)

    # the coefficient is undefined if one or both of the variables are constant => annotate the reason
    if len(x_unique) == 1 or len(y_unique) == 1:
        if x_unique == {0} and y_unique == {0}:
            return np.nan, 'expression and modification constant 0'
        if x_unique == {0}:
            return np.nan, 'expression constant 0'
        if y_unique == {0}:
            return np.nan, 'modification constant 0'

        if len(x_unique) > len(y_unique):
            return np.nan, 'modification constant ({0})'.format(list(y_unique)[0])
        if len(y_unique) > len(x_unique):
            return np.nan, 'expression constant ({0})'.format(list(x_unique)[0])

        assert not (len(x_unique) == 1 and len(y_unique) == 1)

        x_unique = list(x_unique)[0]
        y_unique = list(y_unique)[0]

        if (x_unique, y_unique) in [(1, 1), (-1, -1)]:
            return np.nan, 'expression == modification ({0}, {1})'.format(x_unique, y_unique)

        if (x_unique, y_unique) in [(-1, 1), (1, -1)]:
            return np.nan, 'expression <> modification ({0}, {1})'.format(x_unique, y_unique)

        raise RuntimeError('{0}, {1}'.format(x_unique, y_unique))

    return pearsonr(x, y)[0], 'success'


def is_nan(x: float) -> bool:
    """
    Check if x is a valid number.
    """
    if x is None:
        return True
    return np.isnan(x)


def z_score(r: float) -> float:
    """
    Fisher's correlation coefficient to z-score transform
    :param r: Pearson correlation coefficient
    :return: approximately normal z-score
    """
    return 0.5 * np.log((1 + r) / (1 - r))


def fisher_correlation_comparison(r_1: float, n_1: int, r_2: float, n_2: int) -> float:
    """
    Statistical test for comparing two Pearson correlation coefficients from two independent samples.
    :param r_1: first correlation coefficient
    :param n_1: sample size of the first coefficient
    :param r_2: second coefficient
    :param n_2: sample size of the second coefficient
    :return: two-sided p-value
    """
    # ensure that the given coefficients are numbers
    if is_nan(r_1) or is_nan(r_2):
        return np.nan

    # the sample sizes need to be > 0
    if not n_1 or not n_2:
        return np.nan

    # Fisher's r to z transform for r_1 and r_2
    z_1 = z_score(r_1)
    z_2 = z_score(r_2)

    # standard error
    se = np.sqrt((1 / (n_1 - 3)) + (1 / (n_2 - 3)))

    # z-score of the comparison
    z = (z_1 - z_2) / se

    # transform of z-score to two-sided p-value
    return norm.sf(abs(z)) * 2


def steiger(n: int, r_1: float, r_2: float, r_3: float) -> float:
    """
    Statistical test for comparing two Pearson correlation coefficients from the same sample with one variable in
    common. For example, correlation between X and Y, and correlation between X and Z.
    :param n: sample size
    :param r_1: first correlation coefficient to be compared, e.g. X and Y
    :param r_2: second correlation coefficient to be compared, e.g. X and Z
    :param r_3: correlation coefficient between the second variables, e.g. Y and Z
    :return: two-sided p-value
    """
    # ensure that the given coefficients are numbers
    if is_nan(r_1) or is_nan(r_2) or is_nan(r_3):
        return np.nan

    # the sample size needs to be > 0
    if not n:
        return np.nan

    # Fisher's r to z transform for r_1 and r_2
    z_1 = z_score(r_1)
    z_2 = z_score(r_2)

    # squared mean
    r_ms = pow((r_1 + r_2) * 0.5, 2)

    # error component
    c = (r_3 * (1 - 2 * r_ms) - 0.5 * r_ms * (1 - 2 * r_ms - pow(r_3, 2))) / pow(1 - r_ms, 2)

    # z-score of the comparison
    z = ((z_1 - z_2) * np.sqrt(n - 3)) / np.sqrt(2 - 2 * c)

    # transform of z-score to two-sided p-value
    return norm.sf(abs(z)) * 2
