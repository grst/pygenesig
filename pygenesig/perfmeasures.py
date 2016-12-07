"""
Collection of functions to compute various performance measures from a 2x2 confusion matrix.
As an introduction to evaluating classificators, we recommend reading this `paper`_ about ROC analysis.

.. _paper: http://dx.doi.org/10.1016/j.patrec.2005.10.010.
"""


import numpy as np
from math import sqrt


def mcc(TP, FN, FP, TN):
    """Matthews Correlation Coefficient"""
    return np.divide(TP * TN + FP * FN, sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)))


def sens(TP, FN, FP, TN):
    """Sensitivity"""
    return np.divide(TP, (TP + FN))


def spec(TP, FN, FP, TN):
    """Specificity"""
    return np.divide(TN, (TN + FP))


def prec_pos(TP, FN, FP, TN):
    """Posivitve Precision"""
    return np.divide(TP, (TP + FP))


def recall_pos(TP, FN, FP, TN):
    """Positive recall"""
    return np.divide(TP, (TP + FN))


def prec_neg(TP, FN, FP, TN):
    """Negative Precision"""
    return np.divide(TN, (TN + FN))


def recall_neg(TP, FN, FP, TN):
    """Negative recall"""
    return np.divide(TN, (TN + FP))


def f1_pos(TP, FN, FP, TN):
    """f1-measure on positive instances"""
    return np.divide(2 * prec_pos(TP, FN, FP, TN) * recall_pos(TP, FN, FP, TN),
                     (prec_pos(TP, FN, FP, TN) + recall_pos(TP, FN, FP, TN)))


def f1_neg(TP, FN, FP, TN):
    """f1-measure on negative instances. """
    return np.divide(2 * prec_neg(TP, FN, FP, TN) * recall_neg(TP, FN, FP, TN),
                     (prec_neg(TP, FN, FP, TN) + recall_neg(TP, FN, FP, TN)))
