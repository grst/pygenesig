import unittest
import pandas as pd
import pandas.util.testing as pdt
import numpy as np
from pygenesig.gini import *
from pygenesig.tools import load_gmt, translate_signatures, read_gct
import os


def sort_dict_of_lists(the_dict):
    """Make dict of unordered lists compareable. """
    return {
        key: sorted(value) for key, value in the_dict.items()
    }


class TestGini(unittest.TestCase):
    def setUp(self):
        self.expr = read_gct("./gini_david/data/roche_annotated_exprs.gct")
        self.pdata = pd.read_csv("./gini_david/data/roche_annotated_pdata.tsv", sep="\t")
        self.fdata = pd.read_csv("./gini_david/data/roche_annotated_fdata.tsv", sep="\t")
        self.target = self.pdata.SIG_NAME.as_matrix()

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

        sigs_david_geneid = load_gmt("./gini_david/data/exp.tissuemark.gtex.roche.GeneID.progmt")
        sigs_david_symbol = load_gmt("./gini_david/data/exp.tissuemark.gtex.roche.symbols.gmt")

        sg = GiniSignatureGenerator(self.expr, self.target, min_gini=.8, max_rk=3, min_expr=5)
        sigs = sg.mk_signatures()
        rosetta_geneid = dict(enumerate(self.fdata.GeneID))
        rosetta_symbol = dict(enumerate(self.fdata.GeneSymbol))
        sigs_geneid = translate_signatures(sigs, rosetta_geneid)
        sigs_symbol = translate_signatures(sigs, rosetta_symbol)

        self.assertDictEqual(sort_dict_of_lists(sigs_david_geneid), sort_dict_of_lists(sigs_geneid))
        self.assertDictEqual(sort_dict_of_lists(sigs_david_symbol), sort_dict_of_lists(sigs_symbol))

