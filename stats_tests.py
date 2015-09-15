from scipy.stats import ttest_ind, normaltest


def t_test(series1, series2):
    # Test won't work with NaNs:
    [s.dropna(axis=0, how='any', inplace=True) for s in [series1, series2]]

    t_critical, p_value = ttest_ind(series1, series2)
    return "t-critico = {:.5f}".format(t_critical), \
           "p-value = {:.5f}".format(p_value)


def n_test(series):
    n_test, p_val = normaltest(series)
    return "$p = {:.5f}$".format(p_val)
