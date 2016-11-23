from pygenesig.validation import *
import unittest

class DummySignatureGenerator(SignatureGenerator):
    def mk_signatures(self, subset):
        pass

class DummySignatureTester(SignatureTester):
    def test_signatures(self, subset):
        pass

class TestValidation(unittest.TestCase):
    def test_cross_validation(self):
        pass