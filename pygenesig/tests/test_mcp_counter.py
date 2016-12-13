import unittest
from pygenesig.mcp_counter import *
import logging


class TestMCPCounter(unittest.TestCase):
    def setUp(self):
        self.expr = np.array([2, 3, 5, 4, 9, 15])
        self.target = np.array(["A", "B", "C", "A", "B", "C"])

    def test_fold_change(self):
        self.assertAlmostEqual(-5, fold_change(self.expr, self.target == "A"))

    def test_specific_fold_change(self):
        self.assertAlmostEqual(-3/4, specific_fold_change(self.expr, self.target == "A",
                                                          [self.target == "B", self.target == "C"]))

    def test_roc_auc(self):
        self.assertAlmostEquals(0.875, roc_auc(100-self.expr, self.target == "A"))