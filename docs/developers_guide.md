# Implement own methods 

A major strength of *pygenesig* is its flexibility.
Here, we show how you can implement own methods for generating and testing signatures. 

To standardize the API for creating and testing signatures, pygenesig provides two abstract classes,
`SignatureGenerator` and `SignatureTester` which serve as template for the different methods. 

## Implement a SignatureGenerator
To get started, read the [apidoc](apidoc.html#pygenesig.validation.SignatureGenerator) of 
`SignatureGenerator`. Essentially, all you have to do is to create a child class of 
`SignatureGenerator` and implement the method `_mk_signatures(expr, target)`. 
The method needs to return a signature dictionary as described in the
[data preparation tutorial](prepare_data.html#signature-dictionary).

For example:

```python
class MySignatureGenerator(SignatureGenerator): 
    def _mk_signatures(self, expr, target):
        tissues = set(target)
        return {
            tissue: get_enriched_genes_with_my_method(expr, target, tissue)
                for tissue in tissues
        }
```

You can also override the constructor `__init__()` to pass additional parameters to your class. 
Make always sure to call the parent constructor, though, as it performs some consistency checks on the data: 

```python
class MySignatureGenerator(SignatureGenerator):
    def __init__(self, expr, target, param1=42):
        super(MySignatureGenerator, self).__init__(expr, target)
        self.param1 = param1 

    def _mk_signatures(self, expr, target):
        # ...
```


### Example: GiniSignatureGenerator
This is what our `GiniSignatureGenerator` looks like: 

```python
class GiniSignatureGenerator(SignatureGenerator):
    def __init__(self, expr, target, min_gini=.7, max_rk=3, min_expr=1, max_rel_rk=.33, aggregate_fun=np.median):
        super(GiniSignatureGenerator, self).__init__(expr, target)
        self.min_gini = min_gini
        self.max_rk = max_rk
        self.min_expr = min_expr
        self.aggregate_fun = aggregate_fun
        self.max_rel_rk = max_rel_rk

    def _mk_signatures(self, expr, target):
        df_aggr = collapse_matrix(expr, target, axis=1, aggregate_fun=self.aggregate_fun)
        return get_gini_signatures(df_aggr, min_gini=self.min_gini, max_rk=self.max_rk, min_expr=self.min_expr,
                                   max_rel_rk=self.max_rel_rk)
```

## Implement a SignatureTester
To get started, read the [apidoc](apidoc.html#pygenesig.validation.SignatureTester) of 
`SignatureTester`. Essentially, all you have to do is to create a child class of 
`SignatureTester` and implement the method `_score_signatures(expr, signatures)`.
The method returns a j x n score matrix with j signatures and n samples 

For example, a random predictor, being agnostic of the signatures could look like:
```python
import random 
import numpy as np

class MySignatureTester(SignatureTester):
    def _score_signatures(self, expr, signatures):
        n_tissues = len(signatures) # number of signatures (=j)
        n_samples = expr.shape[1]   # number of samples (=n)
        return np.array ([
            [
                random.random() for i in range(n_samples)
            ] for k in range(n_tissues)
        ])

```

Obviously, you want to implement something meaningful instead!

## Beyond signatures
Actually, instead of passing signature dictionaries, you can pass whatever model between a
`SignatureGenerator` and a `SignatureTester`, as long as the two implementations 
'speak the same language'. For example, one could imagine using a support vector machine
for predicting the tissue. In that case you could pass the trained model instead of the signatures.  
