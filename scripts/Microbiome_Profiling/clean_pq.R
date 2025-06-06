
#' Verify the validity of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Mostly for internal use in MiscMetabar functions.
#'
#' @inheritParams clean_pq
#' @param verbose (logical, default FALSE) If TRUE, prompt some warnings.
#' @param min_nb_seq_sample (numeric) Only used if verbose = TRUE.
#'   Minimum number of sequences per samples to not show warning.
#' @param min_nb_seq_taxa (numeric) Only used if verbose = TRUE.
#'   Minimum number of sequences per taxa to not show warning.
#' @return Nothing if the phyloseq object is valid. An error in the other case.
#'  Warnings if verbose = TRUE
#' @export
#' @author Adrien Taudière
verify_pq <- function(physeq,
                      verbose = FALSE,
                      min_nb_seq_sample = 500,
                      min_nb_seq_taxa = 1) {
  if (!methods::validObject(physeq) ||
      !inherits(physeq, "phyloseq")) {
    stop("The physeq argument is not a valid phyloseq object.")
  }
  if (verbose) {
    if (min(sample_sums(physeq)) < min_nb_seq_sample) {
      warning(
        paste0(
          "At least one of your sample contains less than ",
          min_nb_seq_sample,
          " sequences."
        )
      )
    }
    if (min(sample_sums(physeq)) < min_nb_seq_sample) {
      warning(
        paste0(
          "At least one of your taxa is represent by less than ",
          min_nb_seq_taxa,
          " sequences."
        )
      )
    }
    if (sum(is.na(physeq@sam_data)) > 0) {
      warning("At least one of your samples metadata columns contains NA.")
    }
    if (sum(grepl("^[0-9]", sample_names(physeq)) > 0)) {
      message(
        "At least one sample name start with a number.
      It may introduce bug in some function such
      as psmelt. You may replace sample_names using
      for example :
      sample_names(physeq) <- paste('samp', sample_names(physeq))"
      )
    }
  }
}



################################################################################
#'  Clean phyloseq object by removing empty samples and taxa
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'  In addition, this function check for discrepancy (and rename) between
#' (i) taxa names in refseq, taxonomy table and otu_table and between
#' (ii) sample names in sam_data and otu_table.
#'
#' @param physeq (required): a \code{\link[phyloseq]{phyloseq-class}} object obtained
#'   using the `phyloseq` package.
#' @param remove_empty_samples (logical) Do you want to remove samples
#'   without sequences (this is done after removing empty taxa)
#' @param remove_empty_taxa (logical) Do you want to remove taxa
#'   without sequences (this is done before removing empty samples)
#' @param clean_samples_names (logical) Do you want to clean samples names?
#' @param silent (logical) If true, no message are printing.
#' @param verbose (logical) Additional informations in the message
#'   the verbose parameter overwrite the silent parameter.
#' @param force_taxa_as_columns (logical) If true, if the taxa are rows
#'   transpose the otu_table and set taxa_are_rows to false
#' @param force_taxa_as_rows (logical) If true, if the taxa are columns
#'   transpose the otu_table and set taxa_are_rows to true
#' @param reorder_taxa (logical) if TRUE the otu_table is ordered by the number of
#'   sequences of taxa (ASV, OTU) in descending order. Default to FALSE.
#' @param rename_taxa (logical) if TRUE, taxa (ASV, OTU) are renamed by their position
#'   in the OTU_table and prefix_taxa_names param (by default: Taxa_1, Taxa_2, ...).
#'   Default to FALSE. If rename taxa (ASV, OTU) is true,
#'   the taxa (ASV, OTU) names in verbose information can be misleading.
#' @param simplify_taxo (logical) if TRUE, correct the taxonomy_table using the
#'   `MiscMetabar::simplify_taxo()` function
#' @param prefix_taxa_names (default "Taxa_"): the prefix of taxa names (eg. "ASV_" or "OTU_")
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object
#' @export
#' @author Adrien Taudière
clean_pq <- function(physeq,
                     remove_empty_samples = TRUE,
                     remove_empty_taxa = TRUE,
                     clean_samples_names = TRUE,
                     silent = FALSE,
                     verbose = FALSE,
                     force_taxa_as_columns = FALSE,
                     force_taxa_as_rows = FALSE,
                     reorder_taxa = FALSE,
                     rename_taxa = FALSE,
                     simplify_taxo = FALSE,
                     prefix_taxa_names = "_Taxa") {
  if (clean_samples_names) {
    if (!is.null(physeq@refseq)) {
      if (sum(!names(physeq@refseq) %in% taxa_names(physeq)) > 0) {
        names(physeq@refseq) <- taxa_names(physeq)
        if (!silent) {
          message("Change the samples names in refseq slot")
        }
      }
    }
    if (!is.null(physeq@tax_table)) {
      if (sum(!rownames(physeq@tax_table) %in% taxa_names(physeq)) > 0) {
        rownames(physeq@tax_table) <- taxa_names(physeq)
        if (!silent) {
          message("Change the taxa names in tax_table slot")
        }
      }
    }
    
    if (!is.null(physeq@sam_data)) {
      if (sum(!rownames(physeq@sam_data) %in% sample_names(physeq)) > 0) {
        rownames(physeq@sam_data) <- sample_names(physeq)
        if (!silent) {
          message("Change the samples names in sam_data slot")
        }
      }
    }
  }
  
  verify_pq(physeq)
  
  if (reorder_taxa) {
    physeq <- reorder_taxa_pq(physeq, taxa_names(physeq)[order(taxa_sums(physeq), decreasing = TRUE)])
  }
  
  if (rename_taxa) {
    taxa_names(physeq) <- paste0(prefix_taxa_names, seq(1, ntaxa(physeq)))
  }
  
  if (sum(grepl("^0", sample_names(physeq)) > 0) && !silent) {
    message(
      "At least one sample name start with a zero.
    That can be a problem for some phyloseq functions such as
    plot_bar and psmelt."
    )
  }
  
  if (force_taxa_as_columns && force_taxa_as_rows) {
    stop("You can't force taxa as column and taxa as row in the same time.")
  }
  
  if (force_taxa_as_columns && taxa_are_rows(physeq)) {
    otu_table(physeq) <-
      otu_table(t(as.matrix(unclass(
        physeq@otu_table
      ))), taxa_are_rows = FALSE)
    message("Taxa are now in columns.")
  }
  
  if (force_taxa_as_rows && !taxa_are_rows(physeq)) {
    otu_table(physeq) <-
      otu_table(t(as.matrix(unclass(
        physeq@otu_table
      ))), taxa_are_rows = TRUE)
    message("Taxa are now in rows.")
  }
  
  if (simplify_taxo) {
    physeq <- simplify_taxo(physeq)
  }
  
  new_physeq <- physeq
  
  if (remove_empty_taxa) {
    if (sum(taxa_sums(new_physeq) != 0) > 0) {
      new_physeq <- subset_taxa(physeq, taxa_sums(physeq) > 0)
    }
  }
  if (remove_empty_samples) {
    if (sum(sample_sums(new_physeq) != 0) > 0) {
      new_physeq <- subset_samples(new_physeq, sample_sums(physeq) > 0)
    }
  }
  
  if (verbose) {
    message(
      paste(
        "Cleaning suppress",
        ntaxa(physeq) - ntaxa(new_physeq),
        "taxa (",
        paste(taxa_names(physeq)[taxa_sums(physeq) == 0], collapse = " / "),
        ") and",
        nsamples(physeq) - nsamples(new_physeq),
        "sample(s) (",
        paste(sample_names(physeq)[sample_sums(physeq) == 0], collapse = " / "),
        ")."
      )
    )
  } else if (!silent) {
    message(
      paste(
        "Cleaning suppress",
        ntaxa(physeq) - ntaxa(new_physeq),
        "taxa and",
        nsamples(physeq) - nsamples(new_physeq),
        "samples."
      )
    )
  }
  
  verify_pq(new_physeq)
  return(new_physeq)
}

################################################################################