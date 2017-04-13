import unittest
import numpy as np
import pandas as pd
import random
from pygenesig.tools import *
from pygenesig.file_formats import *
import pandas.util.testing as pdt
from tempfile import NamedTemporaryFile
np.random.seed(42)


def sort_dict_of_lists(the_dict):
    """Make dict of unordered lists compareable. """
    return {
        key: sorted(value) for key, value in the_dict.items()
    }


class TestTools(unittest.TestCase):
    def test_read_write_target(self):
        with NamedTemporaryFile(mode='w') as f:
            data = np.random.choice(["A", "B", "C", "D"], 100)
            write_target(data, f.name)
            read_data = read_target(f.name)
        np.testing.assert_array_equal(data, read_data)

    def test_read_write_rosetta(self):
        with NamedTemporaryFile(mode='w') as f:
            data = np.random.choice(["A", "B", "C", "D"], 100)
            write_rosetta(data, f.name)
            read_data = read_rosetta(f.name, as_dict=False)
        np.testing.assert_array_equal(data, read_data)

    def test_read_write_expr(self):
        with NamedTemporaryFile(mode='w', suffix=".npy") as f:
            data = np.random.rand(2000, 20)
            write_expr(data, f.name)
            read_data = read_expr(f.name)
        np.testing.assert_array_equal(data, read_data)

    def test_read_write_gct(self):
        with NamedTemporaryFile(mode='w') as f:
            data = np.random.rand(2000, 20)
            write_gct(f.name, data)
            read_data = read_gct(f.name)
        np.testing.assert_array_almost_equal(data, read_data)

    def test_read_write_gmt(self):
        with NamedTemporaryFile(mode='w') as f:
            signatures = {
                tissue: np.random.choice(["A", "B", "C", "D", "E", "F"], 3)
                for tissue in ["heart", "liver", "kidney", "brain"]
            }
            write_gmt(signatures, f.name)
            read_signatures = read_gmt(f.name)
        self.assertEqual(sort_dict_of_lists(signatures), sort_dict_of_lists(read_signatures))


