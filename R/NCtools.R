#devtools::install_github("annennenne/causalDisco")
library(causalDisco) 

# Simulate DAG negative control
ncDAG <- function(d, nedges, permute = FALSE, type = "ER") {
  adjm <- matrix(1, d, d)
  adjm[upper.tri(adjm, diag = TRUE)] <- 0
  adjm_lowtri <- which(lower.tri(adjm))
  adjm[sample(adjm_lowtri, length(adjm_lowtri)-nedges)] <- 0
  
  if (permute) {
    perm <- sample.int(d, replace = FALSE)
    adjm <- adjm[perm, perm]
  }
  
  rownames(adjm) <- colnames(adjm) <- paste("x", 1:d, sep = "")
  adjm
}

# Convert DAG adjacency matrix into CPDAG adjacency matrix
as.cpdag <- function(amat) {
  graph2amat(pcalg::dag2cpdag(as.graphNEL(amat)))
}

# Simulate CPDAG negative control
ncCPDAG <- function(d, nedges, permute = TRUE, type = "ER") {
  as.cpdag(ncDAG(d, nedges, permute, type))
}


# Test for overall skeleton fit
# @param conf Should be an adjacency confusion matrix as returned from
# e.g. causalDisco::confusion()
# @value Returns list with the following information: The observed number of 
# true positives, the expected number of true positives under random guessing,
# and a one-sided p-value for their difference. 
skelfit.test <- function(conf) {
  mest <- conf$tp + conf$fp
  mtrue <- conf$tp + conf$fn
  mmax <- conf$tp + conf$fp + conf$fn + conf$tn
  
  p <- phyper(q = conf$tp - 1, 
              m = mtrue,
              n = mmax - mtrue,
              k = mest,
              lower.tail = F)
  list(obs_TP = conf$tp, expected_TP = mest * mtrue / mmax, p = p)
}




# Negative control adjacency precision
nc_adj_precision <- function(mest, mtrue, d, level = 0.95) {
  mmax <- maxnedges(d)
  exp <- mtrue / mmax
  med <- qhyper(0.5, 
                m = mtrue, 
                n = mmax - mtrue, 
                k = mest) / mest
  ci_lwr <- qhyper((1-level)/2, 
                   m = mtrue, 
                   n = mmax - mtrue, 
                   k = mest) / mest
  ci_upr <- qhyper(level + (1-level)/2, 
                   m = mtrue, 
                   n = mmax - mtrue, 
                   k = mest) / mest
  list(expectation = exp, median = med, CI_lwr = ci_lwr, CI_upr = ci_upr)
}


# Negative control adjacency recall
nc_adj_recall <- function(mest, mtrue, d, level = 0.95) {
  mmax <- maxnedges(d)
  exp <- mest / mmax
  med <- qhyper(0.5, 
                m = mtrue, 
                n = mmax - mtrue, 
                k = mest) / mtrue
  ci_lwr <- qhyper((1-level)/2, 
                   m = mtrue, 
                   n = mmax - mtrue, 
                   k = mest) / mtrue
  ci_upr <- qhyper(level + (1-level)/2, 
                   m = mtrue, 
                   n = mmax - mtrue, 
                   k = mest) / mtrue
  list(expectation = exp, median = med, CI_lwr = ci_lwr, CI_upr = ci_upr)
}

# Negative control adjacency F1
nc_adj_f1 <- function(mest, mtrue, d, level = 0.95) {
  mmax <- maxnedges(d)
  exp <- 2 * mest * mtrue / (mmax * mest + mmax * mtrue)
  med <- 2 * qhyper(0.5, 
                    m = mtrue, 
                    n = mmax - mtrue, 
                    k = mest) / (mest + mtrue)
  ci_lwr <- 2 * qhyper((1-level)/2, 
                       m = mtrue, 
                       n = mmax - mtrue, 
                       k = mest) / (mest + mtrue)
  ci_upr <- 2 * qhyper(level + (1-level)/2, 
                       m = mtrue, 
                       n = mmax - mtrue, 
                       k = mest) / (mest + mtrue)
  list(expectation = exp, median = med, CI_lwr = ci_lwr, CI_upr = ci_upr)
}


# Negative control adjacency NPV
nc_adj_npv <- function(mest, mtrue, d, level = 0.95) {
  mmax <- maxnedges(d)
  exp <- 1 - mtrue / mmax
  
  qtransform <- function(q) {
    (mmax - mest - mtrue + q)/(mmax - mest)
  }
  
  med <-  qtransform(qhyper(0.5, 
                            m = mtrue, 
                            n = mmax - mtrue, 
                            k = mest))
  ci_lwr <- qtransform(qhyper((1-level)/2, 
                              m = mtrue, 
                              n = mmax - mtrue, 
                              k = mest))
  ci_upr <- qtransform(qhyper(level + (1-level)/2, 
                              m = mtrue, 
                              n = mmax - mtrue, 
                              k = mest))
  list(expectation = exp, median = med, CI_lwr = ci_lwr, CI_upr = ci_upr)
}

# Negative control adjacency specificity
nc_adj_specificity <- function(mest, mtrue, d, level = 0.95) {
  mmax <- maxnedges(d)
  exp <- 1 - mest / mmax
  
  qtransform <- function(q) {
    (mmax - mest - mtrue + q)/(mmax - mtrue)
  }
  
  med <-  qtransform(qhyper(0.5, 
                            m = mtrue, 
                            n = mmax - mtrue, 
                            k = mest))
  ci_lwr <- qtransform(qhyper((1-level)/2, 
                              m = mtrue, 
                              n = mmax - mtrue, 
                              k = mest))
  ci_upr <- qtransform(qhyper(level + (1-level)/2, 
                              m = mtrue, 
                              n = mmax - mtrue, 
                              k = mest))
  list(expectation = exp, median = med, CI_lwr = ci_lwr, CI_upr = ci_upr)
}

