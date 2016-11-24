import unittest
import pandas as pd
import numpy as np
from pygenesig.gini import *


def sort_dict_of_lists(the_dict):
    """Make dict of unordered lists compareable. """
    return {
        key: sorted(value) for key, value in the_dict.items()
    }


class TestGini(unittest.TestCase):
    def test_aggregate_expression(self):
        expr = np.matrix("0 0 0 0;"
                         "2 4 2 4;"
                         "8 4 6 2")
        target = np.array(["A", "B", "A", "B"])
        expr_aggr = aggregate_expression(expr, target)
        expr_aggr_expected = np.matrix("0 0;"
                                       "2 4;"
                                       "7 3")
        np.testing.assert_array_equal(expr_aggr_expected, expr_aggr)
        subset = np.array([0, 2, 3])
        expr_aggr_subset = aggregate_expression(expr[:, subset], target[subset])
        expr_aggr_subset_expected = np.matrix("0 0;"
                                              "2 4;"
                                              "7 2")
        np.testing.assert_array_equal(expr_aggr_subset_expected, expr_aggr_subset)

    def test_apply_gini_to_dataset(self):
        """Check that the results are consistent with the rogini implementation up to four decimal digits. """
        df_aggr = pd.read_csv("./gini_test.tsv", index_col=0, sep="\t")
        rogini_res = pd.read_csv("./gini_test_0.result.txt", sep="\t")
        gini_res = df_aggr.apply(gini, 1)
        np.testing.assert_almost_equal(list(rogini_res.GINI_IDX), list(gini_res), 4)

    def test_get_gini_signatures(self):
        """Check that the signatures derived with a simple groupby
        from rogini equals the signatures generated by pygenesig.gini. """
        df_aggr = pd.read_csv("./gini_test.tsv", index_col=0, sep="\t")
        rogini_res = pd.read_csv("./gini_test_0-5.result.txt", sep="\t", index_col=0)
        grouped = rogini_res.groupby("CATEGORY")
        rogini_sig = {}
        for key, group in grouped:
            rogini_sig[key] = list(group.index)
        gini_sig = get_gini_signatures(df_aggr, min_gini=.5, max_rk=1, min_expr=0)
        self.assertDictEqual(sort_dict_of_lists(rogini_sig), sort_dict_of_lists(gini_sig))
