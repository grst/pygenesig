import unittest
from pygenesig.bioqc import *
from rpy2.robjects.packages import importr
import logging

base = importr("base")


class TestBioQC(unittest.TestCase):
    def setUp(self):
        # np.matrix does NOT work l
        self.expr = np.array(np.matrix("5. 5. 0. 0.; "
                                       "0. 0. 5. 5.;"
                                       "3. 2. 2. 2.;"
                                       "0. 0. 0. 1."))
        self.target = ["A", "A", "B", "B"]

    def test_signatures2gmt(self):
        signatures = {
            "sig_a": [1, 2, 3, 4],
            "sig_b": [5, 6, 7, 8],
            "sic_c": [9, 10, 11, 12]
        }
        gmt = BioQCSignatureTester.signatures2gmt(signatures)
        self.assertCountEqual(signatures.keys(), list(base.names(gmt)))

    def test_log_pvalue(self):
        """check if log pvalues are identical with running BioQC in R.
        Testcase generated with test_bioqc_log_pvalue.R
        """
        log_p_expected = np.loadtxt("./test_bioqc_log_pvalue.csv", delimiter=',', skiprows=1, usecols=range(1, 5))
        gmt = bioqc.readGmt("./test_bioqc_log_pvalue.gmt")
        gene_symbols = [str(x) for x in range(self.expr.shape[0])]
        p_actual = BioQCSignatureTester.run_bioqc(self.expr, gene_symbols, gmt)
        log_p_actual = -np.log10(p_actual)
        np.testing.assert_array_almost_equal(log_p_expected, log_p_actual)

    def test_test_with_good_signatures(self):
        signatures_good = {
            "A": [0],
            "B": [1]
        }
        tester = BioQCSignatureTester(self.expr, self.target)
        cm = tester.test_signatures(signatures_good, np.array(range(4)))
        cm_expected = np.matrix("2 0;"
                                "0 2") # only TP
        np.testing.assert_array_equal(cm_expected, cm)

    def test_test_with_bad_signatures_and_subset(self):
        signatures_bad = {
            "A": [3, 2],
            "B": [2]
        }
        tester = BioQCSignatureTester(self.expr, self.target)
        cm = tester.test_signatures(signatures_bad, np.array([0, 1, 3]))
        cm_expected = np.matrix("0 2;"
                                "0 1")
        np.testing.assert_array_equal(cm_expected, cm)

    def test_order_of_signatures(self):
        signatures_all = {
            "A": [0],
            "B": [1],
            "C": [3],
            "D": [3],
            "E": [3],
            "F": [3]
        }
        tester = BioQCSignatureTester(self.expr, self.target)
        cm_expected = np.matrix("2 0 0 0 0 0;"
                                "0 2 0 0 0 0;"
                                "0 0 0 0 0 0;"
                                "0 0 0 0 0 0;"
                                "0 0 0 0 0 0;"
                                "0 0 0 0 0 0")
        for i in range(10):
            """test in loop to have random permutations of the dictionary. """
            new_dict = {k: v for k, v in signatures_all.items()}
            cm = tester.test_signatures(new_dict, np.array(range(4)))
            np.testing.assert_array_equal(cm_expected, cm)