from pygenesig.validation import *
import unittest
import tempfile
import numpy as np
import dask.async
from dask.delayed import compute
import dask.threaded
import sklearn.model_selection


class TestValidation(unittest.TestCase):
    class DummySignatureGenerator(SignatureGenerator):
        DUMMY_SIGS = {
            "A": ("A", "B", "C"),
            "B": ("B", "C", "D")
        }

        def _mk_signatures(self, expr, target):
            return self.DUMMY_SIGS

    class DummySignatureTester(SignatureTester):
        PREDICTED = ["A", "B"]

        def _predict(self, expr, signatures):
            return self.PREDICTED

    def setUp(self):
        expr = np.random.random_sample((200, 200))
        self.target = np.array(["A" if x < .3 else "B" for x in np.random.random_sample(200)])
        self.expr_file = tempfile.NamedTemporaryFile()
        self.target_file = tempfile.NamedTemporaryFile()
        np.save(self.expr_file, expr)
        np.savetxt(self.target_file, self.target, delimiter=",", fmt="%s")
        self.expr_file.flush()
        self.target_file.flush()

    def tearDown(self):
        self.expr_file.close()
        self.target_file.close()

    def test_cross_validation(self):
        sig_list, res_list = cross_validate_signatures(self.expr_file.name, self.target_file.name,
                                  self.DummySignatureGenerator, self.DummySignatureTester,
                                  splitter=sklearn.model_selection.StratifiedKFold(n_splits=10))
        # use the non-parallel scheduler for debugging.
        r_sig, r_res = compute(sig_list, res_list, get=dask.async.get_sync)
        self.assertEqual(10, len(r_sig))
        self.assertEqual(10, len(r_res))
        for i in range(10):
            self.assertEqual(self.DummySignatureGenerator.DUMMY_SIGS, r_sig[i])
            actual, predicted = r_res[i]
            # cannot test actual, since it is a subset of self.target
            self.assertEqual(self.DummySignatureTester.PREDICTED, predicted)

    def test_cross_validation_parallel(self):
        sig_list, res_list = cross_validate_signatures(self.expr_file.name, self.target_file.name,
                                                       self.DummySignatureGenerator, self.DummySignatureTester,
                                                       splitter=sklearn.model_selection.StratifiedKFold(n_splits=10))
        # use the threaded scheduler
        r_sig, r_res = compute(sig_list, res_list, get=dask.threaded.get)
        self.assertEqual(10, len(r_sig))
        self.assertEqual(10, len(r_res))
        for i in range(10):
            self.assertEqual(self.DummySignatureGenerator.DUMMY_SIGS, r_sig[i])
            actual, predicted = r_res[i]
            # cannot test actual, since it is a subset of self.target
            self.assertEqual(self.DummySignatureTester.PREDICTED, predicted)
