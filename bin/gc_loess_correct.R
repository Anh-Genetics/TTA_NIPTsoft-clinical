#!/usr/bin/env Rscript
# ============================================================
# gc_loess_correct.R
# EN: GC content and mappability bias correction using LOESS regression.
# VI: Hiệu chỉnh sai lệch GC và tính ánh xạ bằng hồi quy LOESS.
#
# EN: Input : raw bin counts BED, GC content BED, mappability BED
#     Output: normalized counts Parquet, GC bias metrics JSON
# VI: Đầu vào : file BED số đếm bins thô, file BED GC, file BED tính ánh xạ
#     Đầu ra  : file Parquet số đếm chuẩn hóa, file JSON chỉ số sai lệch GC
# ============================================================

suppressPackageStartupMessages({
  library(optparse)   # EN: argument parsing | VI: phân tích đối số dòng lệnh
  library(arrow)      # EN: Parquet I/O     | VI: đọc/ghi file Parquet
  library(jsonlite)   # EN: JSON output     | VI: xuất kết quả JSON
})

# ── Argument parsing / Phân tích đối số ──────────────────────
option_list <- list(
  make_option("--counts",      type="character",
              help="EN: Raw counts BED file | VI: File BED số đếm thô"),
  make_option("--gc-bed",      type="character",
              help="EN: GC content BED file | VI: File BED nội dung GC"),
  make_option("--map-bed",     type="character",
              help="EN: Mappability BED file | VI: File BED tính ánh xạ"),
  make_option("--sample-id",   type="character",
              help="EN: Sample ID | VI: Mã mẫu"),
  make_option("--out-parquet", type="character",
              help="EN: Output Parquet file | VI: File Parquet đầu ra"),
  make_option("--out-metrics", type="character",
              help="EN: Output JSON metrics file | VI: File JSON chỉ số đầu ra"),
  make_option("--min-gc",      type="double",    default=0.20,
              help="EN: Min GC fraction [0.20] | VI: Phân số GC tối thiểu [0.20]"),
  make_option("--max-gc",      type="double",    default=0.80,
              help="EN: Max GC fraction [0.80] | VI: Phân số GC tối đa [0.80]"),
  make_option("--min-map",     type="double",    default=0.50,
              help="EN: Min mappability [0.50] | VI: Tính ánh xạ tối thiểu [0.50]"),
  make_option("--loess-span",  type="double",    default=0.75,
              help="EN: LOESS span [0.75] | VI: Độ trải rộng LOESS [0.75]")
)

opt <- parse_args(OptionParser(option_list=option_list))

message(sprintf("[INFO] EN: Starting GC/mappability correction / VI: Bắt đầu hiệu chỉnh GC/ánh xạ cho mẫu: %s",
                opt[["sample-id"]]))

# ── Load data / Tải dữ liệu ──────────────────────────────────
read_bed <- function(path, col_names) {
  # EN: Read a BED file, ignoring comment lines
  # VI: Đọc file BED, bỏ qua các dòng chú thích
  df <- read.table(path, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(df)[seq_along(col_names)] <- col_names
  df
}

# EN: Read raw counts BED (expected columns: chr, start, end, name, gc, mappability, raw_count)
# VI: Đọc file BED số đếm thô (các cột: chr, start, end, name, gc, mappability, raw_count)
counts_df <- read_bed(opt$counts, c("chr","start","end","name","gc","mappability","raw_count"))

# EN: Merge GC and mappability from separate BED files if provided
# VI: Gộp GC và tính ánh xạ từ các file BED riêng nếu có
if (!is.null(opt[["gc-bed"]])) {
  gc_df  <- read_bed(opt[["gc-bed"]],  c("chr","start","end","gc_content"))
  map_df <- read_bed(opt[["map-bed"]], c("chr","start","end","mappability_val"))
  # EN: Join by coordinate key (chr_start_end)
  # VI: Kết hợp theo khóa tọa độ (chr_start_end)
  key_counts <- paste(counts_df$chr, counts_df$start, counts_df$end, sep="_")
  key_gc     <- paste(gc_df$chr,     gc_df$start,     gc_df$end,     sep="_")
  key_map    <- paste(map_df$chr,    map_df$start,    map_df$end,    sep="_")
  counts_df$gc_frac  <- gc_df$gc_content[match(key_counts, key_gc)]
  counts_df$map_frac <- map_df$mappability_val[match(key_counts, key_map)]
} else {
  # EN: GC and mappability already embedded in the counts BED
  # VI: GC và ánh xạ đã được nhúng trong file BED số đếm
  counts_df$gc_frac  <- counts_df$gc
  counts_df$map_frac <- counts_df$mappability
}

n_total <- nrow(counts_df)

# ── Filter low-quality bins / Lọc bins chất lượng thấp ───────
# EN: Mask bins outside GC/mappability bounds
# VI: Che khuất bins nằm ngoài giới hạn GC/tính ánh xạ
valid_bins <- (
  !is.na(counts_df$gc_frac)  & counts_df$gc_frac  >= opt[["min-gc"]] &
  counts_df$gc_frac  <= opt[["max-gc"]] &
  !is.na(counts_df$map_frac) & counts_df$map_frac >= opt[["min-map"]]
)
counts_df$is_masked <- !valid_bins

# ── GC LOESS correction / Hiệu chỉnh GC bằng LOESS ──────────
# EN: Fit LOESS model on valid bins with non-zero counts, then correct
# VI: Khớp mô hình LOESS trên bins hợp lệ có số đếm > 0, rồi hiệu chỉnh
fit_mask <- valid_bins & counts_df$raw_count > 0

gc_corrected <- rep(NA_real_, nrow(counts_df))

if (sum(fit_mask) > 10) {
  tryCatch({
    message("[INFO] EN: Fitting LOESS model / VI: Đang khớp mô hình LOESS ...")
    loess_fit <- loess(
      raw_count ~ gc_frac,
      data    = counts_df[fit_mask, ],
      span    = opt[["loess-span"]],
      degree  = 2,
      control = loess.control(surface="direct")
    )
    # EN: Predict expected count at each GC value; protect against zero/negative predictions
    # VI: Dự đoán số đếm kỳ vọng tại mỗi giá trị GC; bảo vệ tránh dự đoán 0 hoặc âm
    global_mean  <- mean(counts_df$raw_count[fit_mask], na.rm=TRUE)
    predicted_gc <- predict(loess_fit, newdata=data.frame(gc_frac=counts_df$gc_frac))
    predicted_gc[predicted_gc <= 0 | is.na(predicted_gc)] <- global_mean
    gc_corrected <- counts_df$raw_count * (global_mean / predicted_gc)
    gc_corrected[!valid_bins] <- NA_real_
    message("[INFO] EN: LOESS GC correction applied / VI: Đã áp dụng hiệu chỉnh LOESS GC")
  }, error = function(e) {
    # EN: Fallback to raw counts when LOESS fails
    # VI: Dùng số đếm thô khi LOESS thất bại
    message("WARNING / CẢNH BÁO: LOESS GC correction failed / Hiệu chỉnh LOESS GC thất bại: ",
            conditionMessage(e))
    gc_corrected <<- as.double(counts_df$raw_count)
  })
} else {
  # EN: Too few valid bins for LOESS — use raw counts directly
  # VI: Quá ít bins hợp lệ cho LOESS — dùng số đếm thô trực tiếp
  message("WARNING / CẢNH BÁO: Too few valid bins for LOESS / Quá ít bins hợp lệ cho LOESS; ",
          "using raw counts / dùng số đếm thô")
  gc_corrected <- as.double(counts_df$raw_count)
}

counts_df$gc_corrected_count <- gc_corrected

# ── Mappability correction / Hiệu chỉnh tính ánh xạ ─────────
# EN: Divide GC-corrected counts by mappability score to normalize for
#     regions with lower uniquely-mappable read fractions
# VI: Chia số đếm đã hiệu chỉnh GC cho điểm ánh xạ để chuẩn hóa cho
#     các vùng có phần đọc ánh xạ duy nhất thấp hơn
map_corrected <- gc_corrected / counts_df$map_frac
map_corrected[!valid_bins]            <- NA_real_
map_corrected[is.infinite(map_corrected)] <- NA_real_
counts_df$norm_count <- map_corrected

# ── Within-sample normalization / Chuẩn hóa trong mẫu ───────
# EN: Global median scaling on autosome bins only
# VI: Chia tỷ lệ theo trung vị toàn cục chỉ trên bins NST thường
autosome_bins <- valid_bins & !grepl("^chr[XYM]", counts_df$chr)
global_median <- median(counts_df$norm_count[autosome_bins], na.rm=TRUE)
if (!is.na(global_median) && global_median > 0) {
  counts_df$norm_count <- counts_df$norm_count / global_median
}

# ── QC metrics / Chỉ số QC ───────────────────────────────────
n_valid        <- sum(valid_bins)
n_masked       <- sum(!valid_bins)
# EN: GC residual SD — lower is better (measures remaining GC correlation)
# VI: Độ lệch chuẩn phần dư GC — thấp hơn là tốt hơn (đo tương quan GC còn lại)
gc_residual_sd <- sd(log2(gc_corrected[fit_mask] /
                            mean(gc_corrected[fit_mask], na.rm=TRUE) + 1e-9), na.rm=TRUE)

# EN: Pearson correlation of raw counts vs GC before/after correction
# VI: Tương quan Pearson giữa số đếm thô và GC trước/sau hiệu chỉnh
cor_raw <- if (sum(!is.na(counts_df$gc_frac) & !is.na(counts_df$raw_count)) > 2)
  cor(counts_df$gc_frac[fit_mask], counts_df$raw_count[fit_mask], use="pairwise.complete.obs")
  else NA

cor_corrected <- if (sum(!is.na(counts_df$gc_frac) & !is.na(gc_corrected)) > 2)
  cor(counts_df$gc_frac[fit_mask], gc_corrected[fit_mask], use="pairwise.complete.obs")
  else NA

# ── Write output / Ghi kết quả ───────────────────────────────
# EN: Select output columns and add sample_id
# VI: Chọn cột đầu ra và thêm mã mẫu
out_df <- counts_df[, c("chr","start","end","raw_count","gc_corrected_count",
                         "norm_count","gc_frac","map_frac","is_masked")]
out_df$sample_id <- opt[["sample-id"]]

message("[INFO] EN: Writing normalized counts Parquet / VI: Đang ghi file Parquet số đếm chuẩn hóa ...")
arrow::write_parquet(out_df, opt[["out-parquet"]])

# EN: Write JSON metrics file
# VI: Ghi file JSON chỉ số QC
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

message(sprintf(
  "[DONE / HOÀN THÀNH] EN: GC correction complete / VI: Hoàn thành hiệu chỉnh GC cho mẫu %s: %d valid / %d total bins",
  opt[["sample-id"]], n_valid, n_total
))
