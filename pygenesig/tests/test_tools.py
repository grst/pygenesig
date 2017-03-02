import unittest
import numpy as np
import pandas as pd
from tools import *
import pandas.util.testing as pdt


class TestTools(unittest.TestCase):
    def test_collapse_matrix_0(self):
        mat = np.array(np.matrix("1 1 2 2;"
                                 "4 4 8 8;"
                                 "1 1 2 2;"
                                 "4 4 8 8;"
                                 "1 1 1 1"))
        group_by = ["A", "B", "A", "B", "C"]

        expected = pd.DataFrame(np.matrix("1 1 2 2;"
                                          "4 4 8 8;"
                                          "1 1 1 1"))
        expected.index = ["A", "B", "C"]

        actual = collapse_matrix(mat, group_by, axis=0)
        pdt.assert_frame_equal(actual, expected)

    def test_collapse_matrix_1(self):
        mat = np.array(np.matrix("1 1 2 3;"
                                 "4 4 8 9;"
                                 "1 1 2 3;"
                                 "4 4 8 9;"
                                 "1 1 1 2"))
        group_by = ["A", "B", "A", "B"]

        expected = pd.DataFrame(np.matrix("3 4;"
                                          "12 13;"
                                          "3 4;"
                                          "12 13;"
                                          "2 3"))
        expected.columns = ["A", "B"]

        actual = collapse_matrix(mat, group_by, axis=1, aggregate_fun=np.sum)
        pdt.assert_frame_equal(actual, expected)


