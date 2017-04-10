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

    >>> expr = np.array([2, 3, 5, 4, 9, 15])
    >>> target = np.array(["A", "B", "C", "A", "B", "C"])
    >>> fold_change(expr, target == "A")
    -5.0

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5

    """
    return np.subtract(np.mean(expr[positive_mask]), np.mean(expr[~positive_mask]))


def specific_fold_change(expr, positive_mask, negative_masks):
    """
    Compute the specific fold change of the positive class with respect to all other classes
    for a single gene `G`.

    According to `Becht et al.`_ the specific fold change is defined as

    .. math::
        sFC = (X - \overline{X}_{min})/(\overline{X}_{max} - \overline{X}_{min})

    where :math:`X` is the mean of positive and :math:`\overline{X}_{max}` is the maximum mean over
    all negative classes and :math:`\overline{X}_{min}` is the minimal mean over all negative classes.

    Args:
        expr (np.ndarray): expression of `G` for each sample.
        positive_mask (np.ndarray): boolean mask for `expr` indicating which samples belong to
            the positive class.
        negative_masks (list of np.ndarray): list of boolean masks for `expr` indicating which samples belong
            to the different negative classes.

    Returns:
        float: specific fold change

    >>> expr = np.array([2, 3, 5, 4, 9, 15])
    >>> target = np.array(["A", "B", "C", "A", "B", "C"])
    >>> specific_fold_change(expr, target == 'A', [target == "B", target == "C"])
    -0.75

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
    Compute the Receiver Operator Characteristics Area under the Curve (ROC AUC) for a single gene `G`.
    This tells how well the gene discriminates between the two classes.

    This is a wrapper for the scikit-learn `roc_auc_score`_ function.

    Args:
       expr (np.ndarray): expression of `G` for each sample.
       positive_mask (np.ndarray, dtype=np.bool): boolean mask for `expr` indicating which samples belong to
            the positive class.

    Returns:
        float: roc auc score

    .. _roc_auc_score:
        http://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score

    >>> expr = np.array([2, 3, 5, 4, 9, 15])
    >>> target = np.array(["A", "B", "C", "A", "B", "C"])
    >>> roc_auc(expr, target == "A")
    0.125

    """
    return roc_auc_score(positive_mask, expr)


class MCPSignatureGenerator(SignatureGenerator):
    """
    Implements the procedure described by `Becht et al.`_ for curating signatures.
    A gene is considered a valid `Transcription Marker` (TM) if it meets the following criteria

    * fold change >= ``min_fc`` (default=2)
    * specific fold change >= ``min_sfc`` (default=1.5)
    * AUC ROC >= ``min_auc`` (default=0.97)

    Args:
        expr: m x n gene expression matrix with m genes and n samples.
        target: m-vector with true tissue for each sample
        min_fc: minimal fold change for a gene to be considered as valid TM
        min_sfc: minimal specific fold change for a gene to be considered as valid TM
        min_auc: minimal ROC AUC for a gene to be considered as valid TM

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5

    """

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
        # TODO: this could be sped up significantly by precomputing the means for each class.
        for i in range(expr.shape[0]):
            for cls in classes:
                fc = fold_change(expr[i, :], masks[cls])
                sfc = specific_fold_change(expr[i, :], masks[cls], [mask for k, mask in masks.items() if k != cls])
                auc = roc_auc(expr[i, :], masks[cls])
                if fc >= self.min_fc and sfc >= self.min_sfc and auc >= self.min_auc:
                    signatures[cls].append(i)

        return signatures


class MCPSignatureTester(SignatureTester):
    """
    Implements the `MCPCounter`_ described by Becht et al. in python.

    The principle is super-simple: take the mean of all marker genes as indicator.
    Also see their `R script`_.

    .. _R script:
        https://github.com/ebecht/MCPcounter/blob/a79614eee002c88c64725d69140c7653e7c379b4/Source/R/MCPcounter.R
    """

    def _score_signatures(self, expr, signatures):
        result = np.empty((len(signatures), expr.shape[1]))
        classes = self.sort_signatures(signatures)
        for i, cls in enumerate(classes):
            inds = np.array(signatures[cls])
            for j in range(expr.shape[1]):
                result[i, j] = np.mean(expr[inds, j]) if len(inds) > 0 else np.NAN
        return result
