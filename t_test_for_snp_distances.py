from scipy import stats


def t_test_for_snp_distances(series1, series2):
    # Test won't work with NaNs:
    [s.dropna(axis=0, how='any', inplace=True) for s in [series1, series2]]

    t_critical, p_value = stats.ttest_ind(series1, series2)
    print("'{}' vs. '{}': \n t-crítico = {} \n p-value = {}\n".format(
        series1.name, series2.name, t_critical, p_value))
