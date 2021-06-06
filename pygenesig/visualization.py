"""
This module contains high-level wrapper functions
for visualizing signature validation results.
"""

import pandas as pd
import numpy as np
from pylab import subplots
import seaborn as sns


def heatmap_figsize(nrow, ncol=None):
    """Generate size tuple for heatmap based on number of rows and columns"""
    if ncol is None:
        ncol = nrow
    return ncol * 0.3 + 3, nrow * 0.3


def aggregate_scores(sig_labels, scores, target):
    """
    Aggregate the scores of all samples of the same class.
    This is useful for plotting an average score heatmap.

    Args:
        sig_labels:
        scores:
        target:

    Returns:
        pd.DataFrame: n x m Data Frame with n signatures and m 'mean-of-target-class' columns.

    """
    scores_df = pd.DataFrame(np.transpose(scores))
    scores_df.columns = sig_labels
    scores_df = scores_df.assign(tissue=pd.Series(target))
    return scores_df.groupby("tissue").mean().transpose()


def plot_score_heatmap(heatmap_df, clip=30):
    """
    High-level function to plot a heatmap of a score matrix.

    Args:
        heatmap_df: Dataframe describing the heatmap
        clip: maximum color score (all values above will have the same color)

    Returns:
        fig: figure object
        ax: axis object

    """
    fig, ax = subplots(
        figsize=heatmap_figsize(heatmap_df.shape[0], heatmap_df.shape[1])
    )
    sns.heatmap(
        heatmap_df,
        ax=ax,
        annot=False,
        vmin=0,
        linewidths=0.2,
        vmax=clip,
        cbar_kws={"label": "BioQC score"},
    )
    ax.set_ylabel("signatures")
    ax.set_xlabel("mean of samples")
    return fig, ax


def plot_confusion_matrix_heatmap(confusion_matrix, labels):
    """
    Visualize a confusion matrix as an annotated heatmap.

    Args:
        confusion_matrix: n x n confusion matrix
        labels: n vector containing the axis labels

    Returns:
        fig: figure object
        ax: axis object

    """
    fig, ax = subplots(figsize=heatmap_figsize(confusion_matrix.shape[0]))
    sns.heatmap(
        confusion_matrix,
        ax=ax,
        xticklabels=labels,
        yticklabels=labels,
        annot=True,
        annot_kws={"size": 9},
        fmt=".0f",
        linewidths=0.2,
        vmin=0,
        vmax=100,
        cbar_kws={"label": "percent of class"},
    )
    ax.set_xlabel("predicted")
    ax.set_ylabel("actual")

    return fig, ax
