import unittest
from pygenesig.file_formats import *
import numpy as np
import tempfile


class TestFileFormats(unittest.TestCase):
    def test_gct(self):
        exprs = read_expr(
            "/pstore/data/bioinfo/users/sturmg/gtex-signatures/data_processed/v3/exprs_processed.npy"
        )
        target = read_target(
            "/pstore/data/bioinfo/users/sturmg/gtex-signatures/data_processed/v3/target.csv"
        )
        # read target used on purpose, as we don't want that in a dict, but an array
        rosetta = read_rosetta(
            "/pstore/data/bioinfo/users/sturmg/gtex-signatures/data_processed/v3/rosetta_processed.csv",
            as_dict=False,
        )

        with tempfile.NamedTemporaryFile(delete=False) as tf:
            write_gct(tf.name, exprs, target, rosetta)
            tf.flush()

            print(tf.name)

            exprs_r, target_r, rosetta_r = read_gct(tf.name, include_pfdata=True)

        np.testing.assert_array_almost_equal(exprs, exprs_r)
        np.testing.assert_array_equal(target, target_r)
        np.testing.assert_array_equal(rosetta, rosetta_r)
