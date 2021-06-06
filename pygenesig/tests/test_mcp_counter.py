import unittest
from pygenesig.mcp_counter import *
from pygenesig.validation import *
import logging


class TestMCPCounter(unittest.TestCase):
    def setUp(self):
        # np.matrix does NOT work l
        self.expr = np.array(
            np.matrix("5. 5. 0. 0.; " "0. 0. 5. 5.;" "3. 2. 2. 2.;" "0. 0. 0. 1.")
        )
        self.target = ["A", "A", "B", "B"]
        self.signatures = {  # signatures used in csv and gmt test files.
            "sig_a": [1, 2, 3, 4],
            "sig_b": [5, 6, 7, 8],
            "sig_c": [9, 10, 11, 12],
        }

    def test_roc_auc(self):
        expr = np.array([2, 3, 5, 4, 9, 15])
        target = np.array(["A", "B", "C", "A", "B", "C"])
        self.assertAlmostEquals(0.875, roc_auc(100 - expr, target == "A"))

    def test_score_signatures(self):
        signatures_bad = {"A": [3, 2], "B": [2]}
        st = MCPSignatureTester(self.expr, self.target)
        score = st.score_signatures(signatures_bad)
        score_expected = np.array(np.matrix("1.5 1.  1.  1.5;" "3.  2.  2.  2.  "))
        np.testing.assert_array_almost_equal(score_expected, score)

    def test_test_with_good_signatures(self):
        signatures_good = {"A": [0], "B": [1]}
        tester = MCPSignatureTester(self.expr, self.target)
        actual, predicted = tester._test_signatures(signatures_good, np.array(range(4)))
        np.testing.assert_array_equal(self.target, predicted)

    def test_test_with_bad_signatures_and_subset(self):
        signatures_bad = {"A": [3, 2], "B": [2]}
        tester = MCPSignatureTester(self.expr, self.target)
        actual, predicted = tester._test_signatures(signatures_bad, np.array([0, 1, 3]))
        cm = tester.confusion_matrix(signatures_bad, actual, predicted)
        cm_expected = np.matrix("0 2;" "0 1")
        np.testing.assert_array_equal(cm_expected, cm)

    def test_order_of_signatures(self):
        signatures_all = {"A": [0], "B": [1], "C": [3], "D": [3], "E": [3], "F": [3]}
        tester = MCPSignatureTester(self.expr, self.target)
        cm_expected = np.matrix(
            "2 0 0 0 0 0;"
            "0 2 0 0 0 0;"
            "0 0 0 0 0 0;"
            "0 0 0 0 0 0;"
            "0 0 0 0 0 0;"
            "0 0 0 0 0 0"
        )
        for i in range(10):
            """test in loop to have random permutations of the dictionary. """
            new_dict = {k: v for k, v in signatures_all.items()}
            actual, predicted = tester._test_signatures(new_dict, np.array(range(4)))
            cm = tester.confusion_matrix(signatures_all, actual, predicted)
            np.testing.assert_array_equal(cm_expected, cm)

    def test_empty_signature_all(self):
        signatures_bad = {"A": [], "C": [], "B": []}
        tester = MCPSignatureTester(self.expr, self.target)
        with self.assertRaises(SignatureTesterException):
            actual, predicted = tester._test_signatures(
                signatures_bad, np.array([0, 1, 3])
            )

    def test_empty_signature1(self):
        signatures_bad = {"A": [3], "C": [], "B": []}
        tester = MCPSignatureTester(self.expr, self.target)
        actual, predicted = tester._test_signatures(signatures_bad, np.array([0, 1, 3]))
        cm = tester.confusion_matrix(signatures_bad, actual, predicted)
        cm_expected = np.matrix("2 0 0;" "1 0 0;" "0 0 0")
        np.testing.assert_array_equal(cm_expected, cm)

    def test_empty_signatures_order(self):
        signatures_all = {"A": [], "B": [2], "C": [3], "D": [], "E": [], "F": [3]}
        tester = MCPSignatureTester(self.expr, self.target)
        cm_expected = np.matrix(
            "0 2 0 0 0 0;"
            "0 2 0 0 0 0;"
            "0 0 0 0 0 0;"
            "0 0 0 0 0 0;"
            "0 0 0 0 0 0;"
            "0 0 0 0 0 0"
        )
        for i in range(10):
            """test in loop to have random permutations of the dictionary. """
            new_dict = {k: v for k, v in signatures_all.items()}
            actual, predicted = tester._test_signatures(new_dict, np.array(range(4)))
            cm = tester.confusion_matrix(signatures_all, actual, predicted)
            np.testing.assert_array_equal(cm_expected, cm)

    def test_dtype(self):
        expr = np.array([["1.5", "2.7"], ["2.8", "3.7"]])
        target = np.array(["A", "B"])
        with self.assertRaises(TypeError):
            st = MCPSignatureTester(expr, target)
