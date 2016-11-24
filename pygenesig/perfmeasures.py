import numpy as np
from math import sqrt


def mcc(TP, FN, FP, TN):
    return np.divide(TP * TN + FP * FN, sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)))


def sens(TP, FN, FP, TN):
    return np.divide(TP, (TP + FN))


def spec(TP, FN, FP, TN):
    return np.divide(TN, (TN + FP))


def prec_pos(TP, FN, FP, TN):
    return np.divide(TP, (TP + FP))


def recall_pos(TP, FN, FP, TN):
    return np.divide(TP, (TP + FN))


def prec_neg(TP, FN, FP, TN):
    return np.divide(TN, (TN + FN))


def recall_neg(TP, FN, FP, TN):
    return np.divide(TN, (TN + FP))


def f1_pos(TP, FN, FP, TN):
    return np.divide(2 * prec_pos(TP, FN, FP, TN) * recall_pos(TP, FN, FP, TN),
                     (prec_pos(TP, FN, FP, TN) + recall_pos(TP, FN, FP, TN)))


def f1_neg(TP, FN, FP, TN):
    return np.divide(2 * prec_neg(TP, FN, FP, TN) * recall_neg(TP, FN, FP, TN),
                     (prec_neg(TP, FN, FP, TN) + recall_neg(TP, FN, FP, TN)))
