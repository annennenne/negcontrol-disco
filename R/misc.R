
# Format confidence intervals for Latex tables
ci_latex_info <- function(vals, level = 0.95) {
  paste("$(",paste(format(quantile(vals, 
                                   probs = c((1 - level)/2, level + (1-level)/2)
  ), digits = 2
  ), 
  collapse = ","),
  ")$", sep = "")
}


# Find v-structures
vstructures <- function(amat) {
  d <- nrow(amat)
  out <- list()
  no_combos <- choose(d-1, 2)
  for (j in 1:d) {
    iks <- combn(setdiff(1:d, j), 2)
    for (l in 1:no_combos) {
      i <- iks[1,l]
      k <- iks[2,l]
      if (amat[i,j] == 0 && amat[k,j] == 0) {
        if (amat[j,i] == 1 && amat[j,k] == 1) {
          if (amat[i,k] == 0 && amat[k,i] == 0) {
            out <- c(out, list(c(i, j, k)))
          }
        }
      }
    }
  }
  
  out
}

# Percentage v-structures recovered
percentVstruct <- function(est_amat, true_amat) {
  v_est <- vstructures(est_amat)
  v_est <- c(v_est, lapply(v_est, rev))
  
  v_true <- vstructures(true_amat)
  no_true <- length(v_true)
  v_true <- c(v_true, lapply(v_true, rev))
  
  no_correct <- length(intersect(v_est, v_true))/2
  
  ifelse(no_true > 0, no_correct/no_true, 1)
}
