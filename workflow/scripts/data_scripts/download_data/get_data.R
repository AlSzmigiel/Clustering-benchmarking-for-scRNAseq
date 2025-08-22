args <- commandArgs(trailingOnly = TRUE)

output_file <- args[1]  # First argument: output file path
dataset_name <- args[2]  # Second argument: dataset name

cat("Output File:", output_file, "\n")
cat("Dataset Name:", dataset_name, "\n")

suppressPackageStartupMessages({
    library(scRNAseq)
    library(matrixStats)
    library(AnnotationHub)
    library(BiocFileCache)
    library(dbplyr)
})

library(ExperimentHub)
library(AnnotationHub)


# ---- cache helpers ----
.resolve_hub_cache <- function(hub = c("AnnotationHub", "ExperimentHub")) {
  hub <- match.arg(hub)
  opt_key <- paste0(hub, ".Cache")
  opt <- getOption(opt_key)
  if (!is.null(opt) && nzchar(opt)) return(opt)

  p <- tryCatch(utils::getFromNamespace("R_user_dir", "tools")(hub, "cache"),
                error = function(e) NULL)
  if (!is.null(p)) return(p)

  if (requireNamespace("rappdirs", quietly = TRUE)) {
    return(rappdirs::user_cache_dir(hub))
  }

  file.path(path.expand("~"), ".cache", hub)
}

# delete EVERYTHING in the cache
purge_annotationhub_cache <- function(rebuild = TRUE, verbose = TRUE) {
  cache <- .resolve_hub_cache("AnnotationHub")
  if (verbose) message("[AnnotationHub] cache at: ", cache)

  suppressWarnings({
    if (dir.exists(cache)) {
      if (requireNamespace("BiocFileCache", quietly = TRUE)) {
        bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
        ids <- try(BiocFileCache::bfcquery(bfc, "")$rid, silent = TRUE)
        if (!inherits(ids, "try-error") && length(ids)) {
          BiocFileCache::bfcremove(bfc, rids = ids)
        }
      }
    }
  })

  if (dir.exists(cache)) {
    unlink(cache, recursive = TRUE, force = TRUE)
    if (verbose) message("[AnnotationHub] cache purged.")
  } else if (verbose) {
    message("[AnnotationHub] cache directory not found; nothing to do.")
  }

  # Optionally recreate a fresh cache
  if (isTRUE(rebuild)) {
    if (verbose) message("[AnnotationHub] rebuilding cache with AnnotationHub()…")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    suppressPackageStartupMessages(require(AnnotationHub))
    ah <- AnnotationHub::AnnotationHub()
    if (verbose) message("[AnnotationHub] ready. records: ", length(ah))
    return(invisible(ah))
  }

  invisible(NULL)
}

# purge ExperimentHub the same way 
purge_experimenthub_cache <- function(rebuild = FALSE, verbose = TRUE) {
  cache <- .resolve_hub_cache("ExperimentHub")
  if (verbose) message("[ExperimentHub] cache at: ", cache)

  suppressWarnings({
    if (dir.exists(cache)) {
      if (requireNamespace("BiocFileCache", quietly = TRUE)) {
        bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
        ids <- try(BiocFileCache::bfcquery(bfc, "")$rid, silent = TRUE)
        if (!inherits(ids, "try-error") && length(ids)) {
          BiocFileCache::bfcremove(bfc, rids = ids)
        }
      }
    }
  })

  if (dir.exists(cache)) {
    unlink(cache, recursive = TRUE, force = TRUE)
    if (verbose) message("[ExperimentHub] cache purged.")
  } else if (verbose) {
    message("[ExperimentHub] cache directory not found; nothing to do.")
  }

  if (isTRUE(rebuild)) {
    if (verbose) message("[ExperimentHub] rebuilding cache with ExperimentHub()…")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    suppressPackageStartupMessages(require(ExperimentHub))
    eh <- ExperimentHub::ExperimentHub()
    if (verbose) message("[ExperimentHub] ready. records: ", length(eh))
    return(invisible(eh))
  }

  invisible(NULL)
}


ah <- purge_annotationhub_cache(rebuild = TRUE)
eh <- purge_experimenthub_cache(rebuild = TRUE)


sce <- do.call(dataset_name, list(ensembl = TRUE))

save(sce, file = output_file)
cat("Dataset downloaded and saved to:", output_file, "\n")