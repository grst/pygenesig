from pygenesig.validation import *
import unittest
import tempfile
import numpy as np
from dask.delayed import compute
import dask.threaded
import sklearn.model_selection
import logging


class TestValidation(unittest.TestCase):
    N_GENES = 1000
    N_SAMPELS = 200

    class DummySignatureGenerator(SignatureGenerator):
        DUMMY_SIGS = {"A": ("A", "B", "C"), "B": ("B", "C", "D")}

        def _mk_signatures(self, expr, target):
            return self.DUMMY_SIGS

    class DummySignatureTester(SignatureTester):
        def _score_signatures(self, expr, signatures):
            scores = np.zeros((len(signatures), expr.shape[1]))
            scores[0, :] = 1.0
            return scores

    def setUp(self):
        expr = np.random.random_sample((self.N_GENES, self.N_SAMPELS))
        target = np.array(
            ["A" if x < 0.3 else "B" for x in np.random.random_sample(self.N_SAMPELS)]
        )
        self.expr_file = tempfile.NamedTemporaryFile()
        self.target_file = tempfile.NamedTemporaryFile()
        np.save(self.expr_file, expr)
        np.savetxt(self.target_file, target, delimiter=",", fmt="%s")
        self.expr_file.flush()
        self.target_file.flush()

    def tearDown(self):
        self.expr_file.close()
        self.target_file.close()

    def test_cross_validation(self):
        sig_list, res_list, train_list, test_list = cv_score(
            self.expr_file.name,
            self.target_file.name,
            self.DummySignatureGenerator,
            self.DummySignatureTester,
            splitter=sklearn.model_selection.StratifiedKFold(n_splits=10),
        )
        # use the non-parallel scheduler for debugging.
        r_sig, r_res, r_train, r_test = compute(
            sig_list, res_list, train_list, test_list
        )
        self.assertEqual(10, len(r_sig))
        self.assertEqual(10, len(r_res))
        for i in range(10):
            self.assertEqual(self.DummySignatureGenerator.DUMMY_SIGS, r_sig[i])
            sm = r_res[i]
            self.assertEqual(
                sm.shape, (len(self.DummySignatureGenerator.DUMMY_SIGS), len(r_test[i]))
            )
            self.assertGreater(np.sum(sm), 0)

    def test_cross_validation_parallel(self):
        sig_list, res_list, train_list, test_list = cv_score(
            self.expr_file.name,
            self.target_file.name,
            self.DummySignatureGenerator,
            self.DummySignatureTester,
            splitter=sklearn.model_selection.StratifiedKFold(n_splits=10),
        )
        # use the threaded scheduler
        r_sig, r_res, r_train, r_test = compute(
            sig_list, res_list, train_list, test_list, get=dask.threaded.get
        )
        self.assertEqual(10, len(r_sig))
        self.assertEqual(10, len(r_res))
        for i in range(10):
            self.assertEqual(self.DummySignatureGenerator.DUMMY_SIGS, r_sig[i])
            sm = r_res[i]
            self.assertEqual(
                sm.shape, (len(self.DummySignatureGenerator.DUMMY_SIGS), len(r_test[i]))
            )
            self.assertGreater(np.sum(sm), 0)

    def test_filter_samples(self):
        target = np.array(["A"] * 2 + ["B"] * 10 + ["C"] * 4 + ["D"] * 14)
        exprs = np.random.rand(1000, 30)
        exprs2, target2 = filter_samples(exprs, target, n_splits=10)
        np.testing.assert_array_equal(target2, np.array(["B"] * 10 + ["D"] * 14))
        mask = np.array(list(range(2, 12)) + list(range(16, 30)))
        np.testing.assert_array_equal(exprs2, exprs[:, mask])
