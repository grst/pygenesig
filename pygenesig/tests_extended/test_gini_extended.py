import unittest
import pandas as pd
import pandas.util.testing as pdt
import numpy as np
from pygenesig.gini import *
import os


def sort_dict_of_lists(the_dict):
    """Make dict of unordered lists compareable. """
    return {
        key: sorted(value) for key, value in the_dict.items()
    }


class TestGini(unittest.TestCase):
    def setUp(self):
        self.expr = np.load("./data/gtex_corr_exprs.npy")
        self.target = np.load("./data/target.npy")

    def test_aggratation_chain(self):
        mat_aggr = collapse_matrix(self.expr, self.target, axis=1)
        mat_aggr2 = collapse_matrix(mat_aggr.as_matrix(), mat_aggr.columns, axis=1)
        pdt.assert_frame_equal(mat_aggr, mat_aggr2)

    def test_gini_with_aggegation(self):
        mat_aggr = collapse_matrix(self.expr, self.target, axis=1)
        sg1 = GiniSignatureGenerator(self.expr, self.target)
        sig1 = sort_dict_of_lists(sg1.mk_signatures())
        sg2 = GiniSignatureGenerator(mat_aggr.as_matrix(), mat_aggr.columns)
        sig2 = sort_dict_of_lists(sg2.mk_signatures())
        self.assertDictEqual(sig1, sig2)
