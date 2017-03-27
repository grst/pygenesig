import unittest
import pandas as pd
import pandas.util.testing as pdt
import numpy as np
from pygenesig.gini import *
from pygenesig.tools import load_gmt, translate_signatures, read_gct, jaccard_ind
import os


def sort_dict_of_lists(the_dict):
    """Make dict of unordered lists compareable. """
    return {
        key: sorted(value) for key, value in the_dict.items()
    }

class TestGini(unittest.TestCase):
    def setUp(self):
        self.expr = read_gct("./gini_david/data/roche_annotated_cpm.gct")
        self.pdata = pd.read_csv("./gini_david/data/roche_annotated_pdata.tsv", sep="\t")
        self.fdata = pd.read_csv("./gini_david/data/roche_annotated_fdata.tsv", sep="\t")
        self.target = self.pdata.SIG_NAME.as_matrix()

    def assertSignaturesAlmostEqual(self, signatures1, signatures2, min_jaccard=.95):
        """Check that two signature dictionaries are almost equal.
        The Jaccard index between the two signature dictionaries needs to be
        above a certain threshold"""
        self.assertEquals(signatures1.keys(), signatures2.keys())
        for tissue in signatures1:
            self.assertGreaterEqual(jaccard_ind(signatures1[tissue], signatures2[tissue]), min_jaccard,
                                    msg="Signature '{}' is not almost equal. "
                                        "Sig1 = {}, Sig2 = {}".format(tissue, signatures1[tissue], signatures2[tissue]))

    def test_gini_with_aggegation(self):
        """ it should not make a difference to generate signatures on an already aggregated matrix."""
        mat_aggr = collapse_matrix(self.expr, self.target, axis=1)
        sg1 = GiniSignatureGenerator(self.expr, self.target)
        sig1 = sort_dict_of_lists(sg1.mk_signatures())
        sg2 = GiniSignatureGenerator(mat_aggr.as_matrix(), mat_aggr.columns)
        sig2 = sort_dict_of_lists(sg2.mk_signatures())
        self.assertDictEqual(sig1, sig2)

    def test_gtex_signatures(self):
        """ The signatures generated with this packaage should be identical
        to those generated with David's R script (see folder gini_david). """

        sigs_david_geneid = sort_dict_of_lists(load_gmt("./gini_david/data/exp.tissuemark.gtex.roche.GeneID.progmt"))
        sigs_david_symbol = sort_dict_of_lists(load_gmt("./gini_david/data/exp.tissuemark.gtex.roche.symbols.gmt"))

        sg = GiniSignatureGenerator(self.expr, self.target, min_gini=.8, max_rk=3, min_expr=5)
        sigs = sg.mk_signatures()
        rosetta_geneid = dict(enumerate((str(x) for x in self.fdata.GeneID)))
        rosetta_symbol = dict(enumerate(self.fdata.GeneSymbol))
        sigs_geneid = sort_dict_of_lists(translate_signatures(sigs, rosetta_geneid))
        sigs_symbol = sort_dict_of_lists(translate_signatures(sigs, rosetta_symbol))

        self.assertSignaturesAlmostEqual(sigs_david_geneid, sigs_geneid, min_jaccard=.9)
        self.assertSignaturesAlmostEqual(sigs_david_symbol, sigs_symbol, min_jaccard=.9)

