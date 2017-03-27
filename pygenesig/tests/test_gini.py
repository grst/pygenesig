import unittest
import pandas as pd
import pandas.util.testing as pdt
import numpy as np
from pygenesig.tools import translate_signatures
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
        expr_aggr = collapse_matrix(expr, target, axis=1)
        expr_aggr_expected = np.matrix("0 0;"
                                       "2 4;"
                                       "7 3")
        np.testing.assert_array_equal(expr_aggr_expected, expr_aggr)
        subset = np.array([0, 2, 3])
        expr_aggr_subset = collapse_matrix(expr[:, subset], target[subset], axis=1)
        expr_aggr_subset_expected = np.matrix("0 0;"
                                              "2 4;"
                                              "7 2")
        np.testing.assert_array_equal(expr_aggr_subset_expected, expr_aggr_subset)

    def test_apply_gini_to_dataset(self):
        """Check that the results are consistent with the rogini implementation up to four decimal digits. """
        df_aggr = pd.read_csv("./gini/input/expression_data.tsv", index_col=0, sep="\t")
        rogini_res = pd.read_csv("./gini/output/expression_data_0.0_1.gini", sep="\t")
        gini_res = df_aggr.apply(gini, 1)
        np.testing.assert_almost_equal(list(rogini_res.GINI_IDX), list(gini_res), 4)

    def test_get_gini_signatures(self):
        """Check that the signatures derived with a simple groupby
        from rogini equals the signatures generated by pygenesig.gini. """
        df_aggr = pd.read_csv("./gini/input/expression_data.tsv", index_col=0, sep="\t")
        rogini_res = pd.read_csv("./gini/output/expression_data_0.5_1.gini", sep="\t", index_col=0)
        grouped = rogini_res.groupby("CATEGORY")
        rogini_sig = {}
        for key, group in grouped:
            rogini_sig[key] = list(group.index)
        gini_sig = get_gini_signatures(df_aggr, min_gini=.5, max_rk=1, min_expr=0)
        rosetta = dict(enumerate(df_aggr.index))
        gini_sig = translate_signatures(gini_sig, rosetta)
        self.assertDictEqual(sort_dict_of_lists(rogini_sig), sort_dict_of_lists(gini_sig))

    def test_get_rogini_format(self):
        """Check that the rogini output generated by pygenesig matches the 'original' by Laura. """
        df_aggr = pd.read_csv("./gini/input/expression_data.tsv", index_col=0, sep="\t")
        rogini_res = pd.read_csv("./gini/output/expression_data_0.5_1.gini", sep="\t")
        rogini_res = rogini_res.iloc[:, :5]  # skip last column, as not compulsory.

        gini_format = get_rogini_format(df_aggr, min_gini=.5, max_rk=1, min_expr=0)
        pdt.assert_frame_equal(rogini_res, gini_format, check_dtype=False)

    def test_get_rogini_format2(self):
        """Check that the rogini output generated by pygenesig matches the 'original' by Laura. """
        df_aggr = pd.read_csv("./gini/input/test_rank.tsv", index_col=0, sep="\t")
        rogini_res = pd.read_csv("./gini/output/rank_0.5_3.gini", sep="\t")

        gini_format = get_rogini_format(df_aggr, min_gini=.5, max_rk=3, min_expr=0)
        print(gini_format)
        pdt.assert_almost_equal(rogini_res, gini_format, check_dtype=False)
