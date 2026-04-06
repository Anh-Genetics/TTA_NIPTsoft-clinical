#!/usr/bin/env Rscript
# ============================================================
# gc_loess_correct.R
# GC content and mappability bias correction using LOESS
# Input : raw bin counts BED, GC content BED, mappability BED
# Output: normalized counts Parquet, GC bias metrics JSON
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
})

# ── Argument parsing ─────────────────────────────────────────
option_list <- list(
  make_option("--counts",      type="character", help="Raw counts BED file"),
  make_option("--gc-bed",      type="character", help="GC content BED file"),
  make_option("--map-bed",     type="character", help="Mappability BED file"),
  make_option("--sample-id",   type="character", help="Sample ID"),
  make_option("--out-parquet", type="character", help="Output Parquet file"),
  make_option("--out-metrics", type="character", help="Output JSON metrics file"),
  make_option("--min-gc",      type="double",    default=0.20, help="Min GC fraction [0.20]"),
  make_option("--max-gc",      type="double",    default=0.80, help="Max GC fraction [0.80]"),
  make_option("--min-map",     type="double",    default=0.50, help="Min mappability [0.50]"),
  make_option("--loess-span",  type="double",    default=0.75, help="LOESS span [0.75]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ── Load data ────────────────────────────────────────────────
read_bed <- function(path, col_names) {
  df <- read.table(path, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(df)[seq_along(col_names)] <- col_names
  df
}

counts_df <- read_bed(opt$counts,    c("chr","start","end","name","gc","mappability","raw_count"))
# If GC and mappability are provided separately
if (!is.null(opt[["gc-bed"]])) {
  gc_df  <- read_bed(opt[["gc-bed"]],  c("chr","start","end","gc_content"))
  map_df <- read_bed(opt[["map-bed"]], c("chr","start","end","mappability_val"))
  # Merge by coordinates
  key_counts <- paste(counts_df$chr, counts_df$start, counts_df$end, sep="_")
  key_gc     <- paste(gc_df$chr,     gc_df$start,     gc_df$end,     sep="_")
  key_map    <- paste(map_df$chr,    map_df$start,    map_df$end,    sep="_")
  counts_df$gc_frac  <- gc_df$gc_content[match(key_counts, key_gc)]
  counts_df$map_frac <- map_df$mappability_val[match(key_counts, key_map)]
} else {
  counts_df$gc_frac  <- counts_df$gc
  counts_df$map_frac <- counts_df$mappability
}

n_total <- nrow(counts_df)

# ── Filter low-quality bins ───────────────────────────────────
valid_bins <- (
  !is.na(counts_df$gc_frac)  & counts_df$gc_frac  >= opt[["min-gc"]] &
  counts_df$gc_frac  <= opt[["max-gc"]] &
  !is.na(counts_df$map_frac) & counts_df$map_frac >= opt[["min-map"]]
)
counts_df$is_masked <- !valid_bins

# ── GC LOESS correction ───────────────────────────────────────
# Use valid bins with non-zero counts for model fitting
fit_mask <- valid_bins & counts_df$raw_count > 0

gc_corrected <- rep(NA_real_, nrow(counts_df))

if (sum(fit_mask) > 10) {
  tryCatch({
    loess_fit <- loess(
      raw_count ~ gc_frac,
      data   = counts_df[fit_mask, ],
      span   = opt[["loess-span"]],
      degree = 2,
      control= loess.control(surface="direct")
    )
    # Predict expected count at each GC value
    global_mean <- mean(counts_df$raw_count[fit_mask], na.rm=TRUE)
    predicted_gc <- predict(loess_fit, newdata=data.frame(gc_frac=counts_df$gc_frac))
    predicted_gc[predicted_gc <= 0 | is.na(predicted_gc)] <- global_mean
    gc_corrected <- counts_df$raw_count * (global_mean / predicted_gc)
    gc_corrected[!valid_bins] <- NA_real_
  }, error = function(e) {
    message("WARNING: LOESS GC correction failed: ", conditionMessage(e))
    gc_corrected <<- as.double(counts_df$raw_count)
  })
} else {
  message("WARNING: Too few valid bins for LOESS; using raw counts")
  gc_corrected <- as.double(counts_df$raw_count)
}

counts_df$gc_corrected_count <- gc_corrected

# ── Mappability correction ────────────────────────────────────
map_corrected <- gc_corrected / counts_df$map_frac
map_corrected[!valid_bins] <- NA_real_
map_corrected[is.infinite(map_corrected)] <- NA_real_
counts_df$norm_count <- map_corrected

# ── Within-sample normalization (per-chromosome then global) ──
# Global median scaling
autosome_bins <- valid_bins & !grepl("^chr[XYM]", counts_df$chr)
global_median <- median(counts_df$norm_count[autosome_bins], na.rm=TRUE)
if (!is.na(global_median) && global_median > 0) {
  counts_df$norm_count <- counts_df$norm_count / global_median
}

# ── QC metrics ───────────────────────────────────────────────
n_valid        <- sum(valid_bins)
n_masked       <- sum(!valid_bins)
gc_residual_sd <- sd(log2(gc_corrected[fit_mask] / mean(gc_corrected[fit_mask], na.rm=TRUE) + 1e-9), na.rm=TRUE)

# Correlation of raw count vs GC (lower is better post-correction)
cor_raw <- if(sum(!is.na(counts_df$gc_frac) & !is.na(counts_df$raw_count)) > 2)
  cor(counts_df$gc_frac[fit_mask], counts_df$raw_count[fit_mask], use="pairwise.complete.obs")
  else NA

cor_corrected <- if(sum(!is.na(counts_df$gc_frac) & !is.na(gc_corrected)) > 2)
  cor(counts_df$gc_frac[fit_mask], gc_corrected[fit_mask], use="pairwise.complete.obs")
  else NA

# ── Write output ─────────────────────────────────────────────
out_df <- counts_df[, c("chr","start","end","raw_count","gc_corrected_count",
                         "norm_count","gc_frac","map_frac","is_masked")]
out_df$sample_id <- opt[["sample-id"]]

arrow::write_parquet(out_df, opt[["out-parquet"]])

# JSON metrics
metrics <- list(
  sample_id         = opt[["sample-id"]],
  n_total_bins      = n_total,
  n_valid_bins      = n_valid,
  n_masked_bins     = n_masked,
  gc_cor_raw        = round(cor_raw, 4),
  gc_cor_corrected  = round(cor_corrected, 4),
  gc_residual_sd    = round(gc_residual_sd, 4),
  loess_span        = opt[["loess-span"]],
  global_median     = round(global_median, 6)
)

writeLines(jsonlite::toJSON(metrics, auto_unbox=TRUE, pretty=TRUE), opt[["out-metrics"]])

message(sprintf("GC correction complete for %s: %d valid / %d total bins",
                opt[["sample-id"]], n_valid, n_total))
