import unittest
from pygenesig.bioqc import *
from pygenesig.validation import *
from pygenesig.tools import load_gmt
from rpy2.robjects.packages import importr
import random
from . import TESTDATA

base = importr("base")


class TestBioQC(unittest.TestCase):
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

    def test_signatures2gmt(self):
        gmt = BioQCSignatureTester.signatures2gmt(self.signatures)
        self.assertCountEqual(self.signatures.keys(), list(base.names(gmt)))

    def test_log_pvalue(self):
        """check if log pvalues are identical with running BioQC in R.
        Testcase generated with test_bioqc_log_pvalue.R
        """
        log_p_expected = np.loadtxt(
            TESTDATA / "./bioqc/test_bioqc_log_pvalue.csv",
            delimiter=",",
            skiprows=1,
            usecols=range(1, 5),
        )
        gmt = bioqc.readGmt(str(TESTDATA / "./bioqc/test_bioqc_log_pvalue.gmt"))
        gene_symbols = [str(x) for x in range(self.expr.shape[0])]
        p_actual = BioQCSignatureTester.run_bioqc(self.expr, gene_symbols, gmt)
        log_p_actual = -np.log10(p_actual)
        np.testing.assert_array_almost_equal(log_p_expected, log_p_actual)

    def test_score_signatures(self):
        log_p_expected = np.loadtxt(
            TESTDATA / "./bioqc/test_bioqc_log_pvalue.csv",
            delimiter=",",
            skiprows=1,
            usecols=range(1, 5),
        )
        st = BioQCSignatureTester(self.expr, self.target)
        log_p = st.score_signatures(self.signatures)
        # works, because the GMT is properly sorted.
        np.testing.assert_array_almost_equal(log_p_expected, log_p)

    def test_test_with_good_signatures(self):
        signatures_good = {"A": [0], "B": [1]}
        tester = BioQCSignatureTester(self.expr, self.target)
        actual, predicted = tester._test_signatures(signatures_good, np.array(range(4)))
        np.testing.assert_array_equal(self.target, predicted)

    def test_test_with_bad_signatures_and_subset(self):
        signatures_bad = {"A": [3, 2], "B": [2]}
        tester = BioQCSignatureTester(self.expr, self.target)
        actual, predicted = tester._test_signatures(signatures_bad, np.array([0, 1, 3]))
        cm = tester.confusion_matrix(signatures_bad, actual, predicted)
        cm_expected = np.matrix("0 2;" "0 1")
        np.testing.assert_array_equal(cm_expected, cm)

    def test_order_of_signatures(self):
        signatures_all = {"A": [0], "B": [1], "C": [3], "D": [3], "E": [3], "F": [3]}
        tester = BioQCSignatureTester(self.expr, self.target)
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
        tester = BioQCSignatureTester(self.expr, self.target)
        with self.assertRaises(SignatureTesterException):
            actual, predicted = tester._test_signatures(
                signatures_bad, np.array([0, 1, 3])
            )

    def test_empty_signature1(self):
        signatures_bad = {"A": [3], "C": [], "B": []}
        tester = BioQCSignatureTester(self.expr, self.target)
        actual, predicted = tester._test_signatures(signatures_bad, np.array([0, 1, 3]))
        cm = tester.confusion_matrix(signatures_bad, actual, predicted)
        cm_expected = np.matrix("2 0 0;" "1 0 0;" "0 0 0")
        np.testing.assert_array_equal(cm_expected, cm)

    def test_empty_signatures_order(self):
        signatures_all = {"A": [], "B": [2], "C": [3], "D": [], "E": [], "F": [3]}
        tester = BioQCSignatureTester(self.expr, self.target)
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
            st = BioQCSignatureTester(expr, target)

    def test_reorder(self):
        template = ["A", "B", "C", "D"]
        target = ["D", "B", "C", "A"]
        np.testing.assert_array_equal(
            [3, 1, 2, 0], BioQCSignatureTester._reorder(template, target)
        )

    def test_reorder_with_random_data(self):
        template = np.array(list(range(1000)))
        target = list(range(1000))
        random.shuffle(target)
        target = np.array(target)
        reorder = BioQCSignatureTester._reorder(template, target)
        np.testing.assert_array_equal(template, target[reorder])
