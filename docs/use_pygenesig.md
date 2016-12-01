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

With the actual and predicted labels, we can create a confusion matrix to understand the quality of the prediction. Pygenesig comes with a convenient wrapper method to create the confusion matrix, but essentially you can use whatever performance measure from [scikit-learn](http://scikit-learn.org/stable/modules/classes.html#sklearn-metrics-metrics). 


```python
import seaborn as sns
confmat = st.confusion_matrix(signatures, actual, predicted)
sig_labels = st.sort_signatures(signatures)
sns.heatmap(confmat, xticklabels=sig_labels, yticklabels=sig_labels
```

![heatmap](_static/img/validate_single_heatmap.png)

```eval_rst
.. Note::
    As python dictionaries have no particular order you can use `SignatureTester.sort_signatures` to obtain a reproducable order of the signatures. 
```

# Putting it together: crossvalidation 
To avoid overfitting we can use *crossvalidation* to create and test signatures. The basic concept is to

* create 10 independet, *stratified* folds (*i.e.* every fold contains about the same amount of items from every class)
* run a `SignatureGenerator` on 9 of the 10 folds to generate signatures
* run a `SignatureTester` on the remaining fold
* repeat this such that every fold has been used for testing once
* aggretate the results. 

# Case studies

