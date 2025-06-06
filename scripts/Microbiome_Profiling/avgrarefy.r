# rarefy community matrix by averaging the rarification 100 iterations

# x is a community matrix with taxa as columns (t(otu_table(phyloseq)))

avgrarefy <- function(x, sample, iterations = 100, seed)
{
  set.seed(seed)
  if (is.na(sample))
    stop("invalid subsampling depth")
  if (is.na(iterations))
    stop("invalid iteration count")
  if (sample < min(rowSums(x)))
    stop("rarefaction depth too low; losing samples")
  inputcast <- x
  tablist <- lapply(seq_len(iterations), function(i) {
    # Suppress warnings because it will otherwise return many warnings about
    # subsampling depth not being met, which we deal with below by returning
    # samples that do not meet the threshold.
    inputcast <- suppressWarnings(rrarefy(inputcast, sample = sample))
    # Remove those that did not meet the depth cutoff
    inputcast <- inputcast[c(rowSums(inputcast) %in% sample),]
    as.matrix(inputcast)
  })
  afunc <- array(
    unlist(as.matrix(tablist)),
    c(dim(as.matrix(tablist[[1]])), length(tablist)))
  output <- apply(afunc, 1:2, mean)
  rownames(output) <- rownames(x)
  colnames(output) <- colnames(x)
  output <- round(output, 0)
  return(output)
}
