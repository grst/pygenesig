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

        def mk_signatures(self, subset):
            return self.DUMMY_SIGS

    class DummySignatureTester(SignatureTester):
        DUMMY_CONFUSION_MATRIX = np.matrix("1,2; 3,4")

        def test_signatures(self, signatures, subset):
            return self.DUMMY_CONFUSION_MATRIX

    def setUp(self):
        expr = np.random.random_sample((200, 200))
        target = np.array(["A" if x < .3 else "B" for x in np.random.random_sample(200)])
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
        sig_list, res_list = cross_validate_signatures(self.expr_file.name, self.target_file.name,
                                  self.DummySignatureGenerator, self.DummySignatureTester,
                                  splitter=sklearn.model_selection.StratifiedKFold(n_splits=10))
        # use the non-parallel scheduler for debugging.
        r_sig, r_res = compute(sig_list, res_list, get=dask.async.get_sync)
        self.assertEqual(10, len(r_sig))
        self.assertEqual(10, len(r_res))
        for i in range(10):
            self.assertEqual(self.DummySignatureGenerator.DUMMY_SIGS, r_sig[i])
            self.assertTrue((self.DummySignatureTester.DUMMY_CONFUSION_MATRIX == r_res[i]).all())

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
            self.assertTrue((self.DummySignatureTester.DUMMY_CONFUSION_MATRIX == r_res[i]).all())