import unittest
import numpy as np
import pandas as pd
import random
from pygenesig.tools import *
import pandas.util.testing as pdt


class TestTools(unittest.TestCase):
    def test_collapse_matrix_0(self):
        """test collapse matrix along axis 0"""
        mat = np.array(np.matrix("1 1 2 2;" "4 4 8 8;" "1 1 2 2;" "4 4 8 8;" "1 1 1 1"))
        group_by = ["A", "B", "A", "B", "C"]

        expected = pd.DataFrame(np.matrix("1 1 2 2;" "4 4 8 8;" "1 1 1 1"))
        expected.index = ["A", "B", "C"]

        actual = collapse_matrix(mat, group_by, axis=0)
        pdt.assert_frame_equal(actual, expected, check_dtype=False)

    def test_collapse_matrix_1(self):
        """test collapse matrix along axis 1"""
        mat = np.array(np.matrix("1 1 2 3;" "4 4 8 9;" "1 1 2 3;" "4 4 8 9;" "1 1 1 2"))
        group_by = ["A", "B", "A", "B"]

        expected = pd.DataFrame(np.matrix("3 4;" "12 13;" "3 4;" "12 13;" "2 3"))
        expected.columns = ["A", "B"]

        actual = collapse_matrix(mat, group_by, axis=1, aggregate_fun=np.sum)
        pdt.assert_frame_equal(actual, expected, check_dtype=False)

    def test_collapse_matrix_1_example2(self):
        """test collapse matrix along axis 1 with selecting a subset of columns."""
        expr = np.matrix("0 0 0 0;" "2 4 2 4;" "8 4 6 2")
        target = np.array(["A", "B", "A", "B"])
        expr_aggr = collapse_matrix(expr, target, axis=1)
        expr_aggr_expected = np.matrix("0 0;" "2 4;" "7 3")
        np.testing.assert_array_equal(expr_aggr_expected, expr_aggr)
        subset = np.array([0, 2, 3])
        expr_aggr_subset = collapse_matrix(expr[:, subset], target[subset], axis=1)
        expr_aggr_subset_expected = np.matrix("0 0;" "2 4;" "7 2")
        np.testing.assert_array_equal(expr_aggr_subset_expected, expr_aggr_subset)

    def test_collapse_matrix_chain(self):
        """re-aggregating an already aggregated matrix should return the identity."""
        nrow = 2000
        ncol = 1000
        expr = np.random.rand(nrow, ncol)
        target = np.random.choice(list("ABCDEFGHIJ"), ncol)
        mat_aggr = collapse_matrix(expr, target, axis=1)
        mat_aggr2 = collapse_matrix(mat_aggr.values, mat_aggr.columns, axis=1)
        pdt.assert_frame_equal(mat_aggr, mat_aggr2)
