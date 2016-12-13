"""
`MCPCounter`_ is a method described by Becht et al., Genome Biology (2015)
for deconvolution of tissue-infiltrating immune and stromal cell populations.

Besides presenting a method for gene-expression based deconvolution they put a
lot of effort into curating signatures.

This module

* re-implements MCP counter in python as SignatureTester
* implements a SignatureGenerator based on the MCP method for curating signatures.

.. _MCPCounter:
    http://dx.doi.org/10.1186/s13059-016-1070-5
"""


from pygenesig.validation import SignatureGenerator, SignatureTester
import numpy as np
from sklearn.metrics import roc_auc_score


def fold_change(expr, positive_mask):
    """
    Compute the fold change between the
    positive and negative samples for a single gene `G`.

    According to `Becht et al.`_ the fold change is defined as

    .. math::
        FC = X - \overline{X}

    where :math:`X` is the mean of positive and :math:`\overline{X}` is the mean of the negative samples.

    Args:
        expr (np.ndarray): expression of `G` for each sample.
        positive_mask (np.ndarray, dtype=np.bool): boolean mask for `expr` indicating which samples belong to
            the positive class.

    Returns:
        float: fold change

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5

    """
    return np.subtract(np.mean(expr[positive_mask]), np.mean(expr[~positive_mask]))


def specific_fold_change(expr, positive_mask, negative_masks):
    """
    Compute the specific fold change of the positive class with respect to all other classes.

    According to `Becht et al.`_ the specific fold change is defined as

    .. math::
        sFC = X - \overline{X}

    where :math:`X` is the mean of positive and :math:`\overline{X}` is the mean of the negative samples.

    Args:
        expr (np.ndarray):
        positive_mask (np.ndarray):
        negative_masks (list of np.ndarray):

    Returns:
        float: specific fold change

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5

    """
    mean_per_class = [np.mean(expr[class_inds]) for class_inds in negative_masks]
    x_min = np.min(mean_per_class)
    x_max = np.max(mean_per_class)
    return np.divide(
        np.subtract(np.mean(expr[positive_mask]), x_min),
        np.subtract(x_max, x_min)
    )


def roc_auc(expr, positive_mask):
    """

    Args:
       expr (np.ndarray):
       positive_mask (np.ndarray, dtype=np.bool):

    Returns:

    """
    return roc_auc_score(positive_mask, normalize(expr))


class MCPSignatureGenerator(SignatureGenerator):
    def __init__(self, expr, target, min_fc=2, min_sfc=1.5, min_auc=.97):
        super(MCPSignatureGenerator, self).__init__(expr, target)
        self.min_fc = min_fc
        self.min_sfc = min_sfc
        self.min_auc = min_auc

    def _mk_signatures(self, expr, target):
        classes = list(set(target))
        masks = {
            cls: target == cls for cls in classes
        }
        signatures = {
            cls: [] for cls in classes
        }
        for i in range(expr.shape[0]):
            for cls in classes:
                fc = fold_change(expr[i, :], masks[cls])
                sfc = specific_fold_change(expr[i, :], masks[cls], [mask for k, mask in masks.items() if k != cls])
                auc = roc_auc(expr[i, :], masks[cls])
                if fc >= self.min_fc and sfc >= self.min_sfc and auc >= self.min_auc:
                    signatures[cls].append(i)

        return signatures


class MCPSignatureTester(SignatureTester):
    def _predict(self, expr, signatures):
        predicted = []
        classes = list(iter(signatures.keys()))
        for j in range(expr.shape[1]):
            sample_means = []
            for cls in classes:
                inds = np.array(signatures[cls])
                cls_mean = np.mean(expr[inds, j])
                sample_means.append(cls_mean)
            predicted.append(classes[np.argmax(sample_means)])
        return predicted


