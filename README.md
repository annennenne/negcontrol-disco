# Negative controls for causal discovery evaluation

This repository contains code for the article

Petersen, A. H. (2025). Are you doing better than random guessing? A call for using negative controls when evaluating causal discovery algorithms. Accepted for *Uncertainty in Artificial Intelligence* 2025. 

A preprint is available on [arXiv](https://arxiv.org/abs/2412.10039).

The file "R/NCtools.R" contains general tools for computing negative controls for causal discovery evaluation.

## Negative control DAGs and CPDAGs. 
 A negative control DAG (Erdős-Rényi type) over *d = 5* nodes with *mest = 7* edges is simulated like this: 

```
source("R/NCtools.R")
ncDAG(d = 5, nedges = 7)
```

and a CPDAG is obtained from

```
source("R/NCtools.R")
ncCPDAG(d = 5, nedges = 7)
```

## Expected negative control values for adjacency metrics 
For adjacency precision, recall, F1 score, negative predictive value (NPV) and specificity, we provide easy computation of negative control/random guessing expected values, medians and confidence intervals. For example, for an estimated graph with *d = 5* nodes, true number of edges *mtrue = 8* and estimated number of edges *mest = 7*, we find the negative control adjacency precision expectation, median and 95% confidence interval like this:

```
nc_adj_recall(mest = 7, mtrue = 8, d = 5, level = 0.95)
```

Similar functions with the same syntax exist for other adjacency metrics; `nc_adj_precision()`, `nc_adj_f1()`, `nc_adj_npv()` and `nc_adj_specificity()`.

## Test for overall skeleton fit
For a given adjacency confusion matrix (with TP = 6, FP = 1, FN = 2 and TN = 1), we conduct a test of overall skeleton fit as follows:

```
# Construct confusion matrix
thisconfusion <- list(tp = 6, fp = 1, fn = 2, tn = 1)

# Test of overall skeleton fit
skelfit.test(thisconfusion)
```
The test function returns a one-sided p-value. 

