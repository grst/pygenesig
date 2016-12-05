# Creating signatures
*pygenesig* shippes with two classes for signature generation. The GiniSignatureGenerator and the LimmaSignatureGenerator. The two methods are described...

First, we load the expression data and target annotation we [prepared earlier](prepare_data.html):
```python
import numpy as np
expr = np.load("expression_matrix.npy")
target = np.genfromtxt("target.csv", delimiter=",", dtype=str)
```

Now, we can generate signatures with the signature generator of our choice:
```python
sg = GiniSignatureGenerator(expr, target)
signatures = sg.mk_signatures()
```

Which will result in something like
```
{'Adipose Tissue': {81,
  82,
  250,
  304,
  309,
  ...}
```

# Testing signatures
*pygenesig* shippes with the `BioQCSignatureTester`. The method is described in more detail [here].

To test signatures, we initalize the tester with the gene expression data and the target labels. Then, we can test different signature sets on the data:

```python
st = BioQCSignatureTester(expr, target)
actual, predicted = st.test_signatures(signatures)
```

From the list of actual and predicted labels, we can for example create a confusion matrix. *Pygenesig* provides a convenient wrapper method to create the confusion matrix, but essentially you can use whatever performance measure from [scikit-learn](http://scikit-learn.org/stable/modules/classes.html#sklearn-metrics-metrics). 


```python
import seaborn as sns
confmat = st.confusion_matrix(signatures, actual, predicted)
sig_labels = st.sort_signatures(signatures)
sns.heatmap(confmat, xticklabels=sig_labels, yticklabels=sig_labels
```

![heatmap](_static/img/validate_single_heatmap.png)

```eval_rst
.. Note::
    As python dictionaries have no particular order you can use ``SignatureTester.sort_signatures()`` to obtain a reproducable order of the signatures. 
```

# Putting it together: crossvalidation 
To avoid overfitting sample specific noise, we can use *crossvalidation* to create and test signatures. To this end, we divide our data into 10 independent, *stratified* folds (*i.e.* every fold contains about the same amount of items from every class). We always use 9 of the 10 folds for generating the signatures and apply them to the remaining fold for testing. This procedure is illustrated in the following flowchart:

<!-- edit flowchart on https://www.draw.io/?chrome=0&lightbox=1&edit=https%3A%2F%2Fwww.draw.io%2F%23G0BxECzhdeMGwJQXB5ZjNHckRWRzQ&nav=1#G0BxECzhdeMGwJQXB5ZjNHckRWRzQ --> 

![flowchart](_static/img/pygenesig_xval.svg)

# Case studies
We have performed several case studies using *pygenesig* on the [GTEx](http://www.gtexportal.org/home/) dataset. These studies can be understood as 'extended examples' of how to use *pygenesig* and are available on [github](https://github.com/grst/gene-set-study/tree/master/notebooks). 

* [Signature generation and cross-validation](https://github.com/grst/gene-set-study/blob/master/notebooks/validate_gini.ipynb)
* [Grid search for parameter optimization](https://github.com/grst/gene-set-study/blob/master/notebooks/gini-gridsearch.ipynb): We systematically tested different values for the gini parameters `min_gini` and `max_rk`. We found, that gini-index is a robust method for signature generation over a wide range of parameters. 
* [Cross-platform and cross-species validation](https://github.com/grst/gene-set-study/blob/master/notebooks/validate-mouse.ipynb): In order to demonstrate, that the gini-method is robust over different organisms and platforms, we generated gene signatures on the GTEx dataset (human, next generation sequenceing) and applied them to a mouse dataset (Affymetrics microarray). For most of the tissues, the signatures are still able to identify their respective tissue.  

