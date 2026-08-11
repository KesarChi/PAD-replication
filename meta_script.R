#!/usr/bin/env Rscript

###
# @Author: Wanqing Yin && wanqingyin7@163.com
# @Date: 2026-6-13 16:24
# @LastEditors: Wanqing Yin && wanqingyin7@163.com
# @LastEditTime: 2026-6-13 16:24
# @FilePath: ~/Desktop/PAD/Scripts/
# @Description: Institute: Beijing Institute Of Genomics, CAS.
#
# Copyright (c) 2026 by ${wanqingyin7@163.com}, All Rights Reserved.
###

options(stringsAsFactors = FALSE)

usage <- function() {
  cat(
"Usage:
  Rscript meta_script.R -i file1,file2[,file3...] -N n_cols --or-ci orci_cols -p p_cols [--rows rows] [--keep keep_cols] [-o out.tsv]
  Rscript meta_script.R -i file1,file2[,file3...] --case case_cols --ctrl ctrl_cols -b beta_cols -se se_cols -p p_cols [--rows rows] [--keep keep_cols] [-o out.tsv]

Required:
  -i, --input     Comma-separated input files. Files may be CSV or TSV.

Sample size input, choose one mode:
  -N, --n         Comma-separated N column numbers, one per input file.
  --case          Comma-separated case column numbers, one per input file.
  --ctrl          Comma-separated control column numbers, one per input file.

Effect input, choose one mode:
  --or-ci, --or   Comma-separated OR(95%CI) column numbers, one per input file.
  -b, --beta      Comma-separated beta/log(OR) column numbers, one per input file.
  -se, --se       Comma-separated SE column numbers, one per input file.

Optional:
  -p, --p         Comma-separated p-value column numbers, one per input file.
  --rows, -r      DATA row numbers after the header.
                  Multi-file mode: one row per input file, e.g. 1,2,5 means file1 row 1 + file2 row 2 + file3 row 5.
                  Single-file mode: multiple rows from that one file are meta-analyzed, e.g. -i file.tsv -r 1,2,5.
                  Default: 1st data row.
  --keep          Columns from the FIRST input file to keep in output even if not used for meta.
                  Example: --keep 1 keeps the first column from file1. Multiple columns: --keep 1,2.
                  Default: keep column 1 from file1.
  --sep           Delimiter: auto, csv, tsv, comma, tab, or a literal delimiter. Default: auto.
  -o, --out       Output TSV. Default: pheno_meta_result.tsv.
  --study-labels  Comma-separated study labels. Default: file basenames.
  --full          Output all debugging/source columns. Default output is concise.
  -h, --help      Show this help.

Examples:
  Rscript meta_script.R -i a.tsv,b.csv,c.tsv -N 2,2,2 --or-ci 3,3,3 -p 4,4,4 --rows 1,2,5 --keep 1
  Rscript meta_script.R -i a.tsv,b.tsv --case 2,2 --ctrl 3,3 -b 4,4 -se 5,5 -p 6,6 -r 3,7 --keep 1,2
", sep = "")
}

stop_usage <- function(msg) {
  cat("Error: ", msg, "\n\n", sep = "", file = stderr())
  usage()
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
full_output <- any(args %in% c("--full"))
args <- args[!args %in% c("--full")]
if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = 0)
}

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^-", key)) stop_usage(paste0("unexpected positional argument: ", key))
    if (i == length(args)) stop_usage(paste0("missing value after ", key))
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

opts <- parse_args(args)

get_opt <- function(keys, default = NULL) {
  for (key in keys) {
    if (!is.null(opts[[key]])) return(opts[[key]])
  }
  default
}

split_csv <- function(x, name, allow_empty = FALSE) {
  if (is.null(x) || !nzchar(trimws(x))) {
    if (allow_empty) return(character())
    stop_usage(paste0("missing ", name))
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_int_list <- function(x, name) {
  values <- split_csv(x, name)
  out <- suppressWarnings(as.integer(values))
  if (any(is.na(out)) || any(out < 1)) {
    stop_usage(paste0(name, " must contain positive integer column numbers"))
  }
  out
}

require_len <- function(x, n_files, name) {
  if (length(x) != n_files) {
    stop_usage(sprintf("%s length (%d) must equal number of input files (%d)", name, length(x), n_files))
  }
  x
}

parse_rows <- function(x, n_files) {
  if (is.null(x) || !nzchar(trimws(x))) return(rep(1L, n_files))
  rows <- parse_int_list(x, "--rows")
  require_len(rows, n_files, "--rows")
}

parse_keep_cols <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(1L)
  cols <- parse_int_list(x, "--keep")
  unique(cols)
}

normalize_sep <- function(x) {
  x <- if (is.null(x)) "auto" else x
  switch(tolower(x),
    "auto" = "auto",
    "csv" = ",",
    "comma" = ",",
    "tsv" = "\t",
    "tab" = "\t",
    x
  )
}

read_table_auto <- function(path, sep_mode) {
  if (!file.exists(path)) stop("input file does not exist: ", path, call. = FALSE)
  sep <- sep_mode
  if (identical(sep, "auto")) {
    ext <- tolower(tools::file_ext(path))
    sep <- if (ext == "csv") "," else "\t"
  }
  read.table(
    path,
    sep = sep,
    header = TRUE,
    colClasses = "character",
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    fill = TRUE,
    na.strings = c("", "NA", "NaN", "null", "NULL")
  )
}

infer_p_col <- function(dat) {
  cleaned <- tolower(gsub("[^a-z0-9]+", "", names(dat)))
  candidates <- c("pvalue", "pval", "p")
  hit <- match(candidates, cleaned, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) NA_integer_ else hit[[1]]
}

extract_col <- function(dat, col, path) {
  if (col > ncol(dat)) {
    stop(sprintf("file %s has %d columns; requested column %d", path, ncol(dat), col), call. = FALSE)
  }
  dat[[col]]
}

to_num_scalar <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1 || is.na(value) || !is.finite(value)) {
    stop(label, " is not a finite numeric value: ", as.character(x), call. = FALSE)
  }
  value
}

parse_or_ci <- function(x) {
  raw <- trimws(as.character(x))
  raw <- chartr("（）", "()", raw)
  raw <- gsub("\\s+", "", raw)
  raw <- gsub("\u2013|\u2014|~|to", "-", raw, ignore.case = TRUE)
  raw <- gsub(",", "", raw)
  pattern <- "^([0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\(([0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)-([0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\)$"
  m <- regexec(pattern, raw, perl = TRUE)
  z <- regmatches(raw, m)[[1]]
  if (length(z) != 4) stop("cannot parse OR(95%CI): ", raw, call. = FALSE)
  vals <- as.numeric(z[2:4])
  names(vals) <- c("OR", "CI_lower", "CI_upper")
  if (any(!is.finite(vals)) || any(vals <= 0) || vals[["CI_lower"]] > vals[["CI_upper"]]) {
    stop("invalid OR(95%CI): ", raw, call. = FALSE)
  }
  vals
}

format_meta_p_value <- function(x) {
  ifelse(is.na(x), NA_character_, as.character(signif(x, 8)))
}

files <- split_csv(get_opt(c("-i", "--input")), "-i/--input")
n_files <- length(files)

n_cols_raw <- get_opt(c("-N", "--n"))
case_cols_raw <- get_opt(c("--case", "--cases"))
ctrl_cols_raw <- get_opt(c("--ctrl", "--control", "--controls"))
has_n_mode <- !is.null(n_cols_raw)
has_case_ctrl_mode <- !is.null(case_cols_raw) || !is.null(ctrl_cols_raw)
if (has_n_mode && has_case_ctrl_mode) stop_usage("use either -N/--n or --case/--ctrl, not both")
if (!has_n_mode && !has_case_ctrl_mode) stop_usage("provide either -N/--n or both --case and --ctrl")
if (has_case_ctrl_mode && (is.null(case_cols_raw) || is.null(ctrl_cols_raw))) {
  stop_usage("case/control mode requires both --case and --ctrl")
}

n_cols <- if (has_n_mode) require_len(parse_int_list(n_cols_raw, "-N/--n"), n_files, "-N/--n") else rep(NA_integer_, n_files)
case_cols <- if (has_case_ctrl_mode) require_len(parse_int_list(case_cols_raw, "--case"), n_files, "--case") else rep(NA_integer_, n_files)
ctrl_cols <- if (has_case_ctrl_mode) require_len(parse_int_list(ctrl_cols_raw, "--ctrl"), n_files, "--ctrl") else rep(NA_integer_, n_files)

p_cols_raw <- get_opt(c("-p", "--p", "--p-value", "--pvalue"))
p_cols <- if (is.null(p_cols_raw)) rep(NA_integer_, n_files) else require_len(parse_int_list(p_cols_raw, "-p/--p"), n_files, "-p/--p")

or_cols_raw <- get_opt(c("--or-ci", "--or_ci", "--or", "--orci"))
beta_cols_raw <- get_opt(c("-b", "--beta"))
se_cols_raw <- get_opt(c("-se", "--se"))
has_or_mode <- !is.null(or_cols_raw)
has_beta_mode <- !is.null(beta_cols_raw) || !is.null(se_cols_raw)
if (has_or_mode && has_beta_mode) stop_usage("use either OR(95%CI) mode or beta/se mode, not both")
if (!has_or_mode && !has_beta_mode) stop_usage("provide either --or-ci columns or both -b and -se columns")

or_cols <- if (has_or_mode) require_len(parse_int_list(or_cols_raw, "--or-ci"), n_files, "--or-ci") else rep(NA_integer_, n_files)
beta_cols <- if (!is.null(beta_cols_raw)) require_len(parse_int_list(beta_cols_raw, "-b/--beta"), n_files, "-b/--beta") else rep(NA_integer_, n_files)
se_cols <- if (!is.null(se_cols_raw)) require_len(parse_int_list(se_cols_raw, "-se/--se"), n_files, "-se/--se") else rep(NA_integer_, n_files)
if (has_beta_mode && (is.null(beta_cols_raw) || is.null(se_cols_raw))) {
  stop_usage("beta/se mode requires both -b/--beta and -se/--se")
}

raw_rows <- get_opt(c("--rows", "-r"))
if (n_files == 1) {
  rows <- if (is.null(raw_rows) || !nzchar(trimws(raw_rows))) 1L else parse_int_list(raw_rows, "--rows")
} else {
  rows <- parse_rows(raw_rows, n_files)
}
keep_cols <- parse_keep_cols(get_opt(c("--keep")))
sep_mode <- normalize_sep(get_opt(c("--sep"), "auto"))
out_file <- get_opt(c("-o", "--out"), "pheno_meta_result.tsv")

tables <- lapply(files, read_table_auto, sep_mode = sep_mode)
if (is.null(p_cols_raw)) {
  p_cols <- vapply(tables, infer_p_col, integer(1))
}
if (n_files == 1) {
  unit_file_idx <- rep(1L, length(rows))
  unit_rows <- rows
  default_labels <- paste0(basename(files[[1]]), "_row", unit_rows)
} else {
  unit_file_idx <- seq_along(files)
  unit_rows <- rows
  default_labels <- basename(files)
}
n_units <- length(unit_rows)
study_labels <- split_csv(get_opt(c("--study-labels"), paste(default_labels, collapse = ",")), "--study-labels")
study_labels <- require_len(study_labels, n_units, "--study-labels")

expand_cols <- function(cols, name) {
  if (n_files == 1) rep(cols[[1]], n_units) else cols
}
n_cols_unit <- expand_cols(n_cols, "-N/--n")
case_cols_unit <- expand_cols(case_cols, "--case")
ctrl_cols_unit <- expand_cols(ctrl_cols, "--ctrl")
p_cols_unit <- expand_cols(p_cols, "-p/--p")
or_cols_unit <- expand_cols(or_cols, "--or-ci")
beta_cols_unit <- expand_cols(beta_cols, "-b/--beta")
se_cols_unit <- expand_cols(se_cols, "-se/--se")

study_dat <- vector("list", n_units)
meta_row_id <- 1L
row_spec <- paste(unit_rows, collapse = ",")

for (unit_idx in seq_len(n_units)) {
  file_idx <- unit_file_idx[[unit_idx]]
  dat <- tables[[file_idx]]
  path <- files[[file_idx]]
  row_idx <- unit_rows[[unit_idx]]
  if (row_idx > nrow(dat)) {
    stop(sprintf("file %s has %d rows; requested row %d", path, nrow(dat), row_idx), call. = FALSE)
  }

  if (has_n_mode) {
    case_value <- NA_real_
    ctrl_value <- NA_real_
    n_value <- to_num_scalar(extract_col(dat, n_cols_unit[[unit_idx]], path)[[row_idx]], sprintf("%s row %d N", path, row_idx))
  } else {
    case_value <- to_num_scalar(extract_col(dat, case_cols_unit[[unit_idx]], path)[[row_idx]], sprintf("%s row %d case", path, row_idx))
    ctrl_value <- to_num_scalar(extract_col(dat, ctrl_cols_unit[[unit_idx]], path)[[row_idx]], sprintf("%s row %d ctrl", path, row_idx))
    if (case_value < 0 || ctrl_value < 0) stop(sprintf("%s row %d case/ctrl must be >= 0", path, row_idx), call. = FALSE)
    n_value <- case_value + ctrl_value
  }
  if (n_value <= 0) stop(sprintf("%s row %d N must be > 0", path, row_idx), call. = FALSE)

  p_raw <- if (is.na(p_cols_unit[[unit_idx]])) {
    NA_character_
  } else {
    as.character(extract_col(dat, p_cols_unit[[unit_idx]], path)[[row_idx]])
  }
  p_value <- if (is.na(p_cols_unit[[unit_idx]])) {
    NA_real_
  } else {
    p_tmp <- to_num_scalar(p_raw, sprintf("%s row %d p", path, row_idx))
    if (p_tmp < 0 || p_tmp > 1) stop(sprintf("%s row %d p must be between 0 and 1", path, row_idx), call. = FALSE)
    p_tmp
  }

  if (has_or_mode) {
    orci_raw <- extract_col(dat, or_cols_unit[[unit_idx]], path)[[row_idx]]
    orci <- parse_or_ci(orci_raw)
    beta <- log(orci[["OR"]])
    se <- (log(orci[["CI_upper"]]) - log(orci[["CI_lower"]])) / (2 * qnorm(0.975))
    effect_text <- as.character(orci_raw)
    or_value <- orci[["OR"]]
    ci_lower <- orci[["CI_lower"]]
    ci_upper <- orci[["CI_upper"]]
  } else {
    beta <- to_num_scalar(extract_col(dat, beta_cols_unit[[unit_idx]], path)[[row_idx]], sprintf("%s row %d beta", path, row_idx))
    se <- to_num_scalar(extract_col(dat, se_cols_unit[[unit_idx]], path)[[row_idx]], sprintf("%s row %d se", path, row_idx))
    if (se <= 0) stop(sprintf("%s row %d se must be > 0", path, row_idx), call. = FALSE)
    or_value <- exp(beta)
    ci_lower <- exp(beta - qnorm(0.975) * se)
    ci_upper <- exp(beta + qnorm(0.975) * se)
    effect_text <- sprintf("%.2f(%.2f-%.2f)", or_value, ci_lower, ci_upper)
  }

  study_dat[[unit_idx]] <- data.frame(
    meta_row = meta_row_id,
    file_index = file_idx,
    study = study_labels[[unit_idx]],
    file = path,
    source_row = row_idx,
    source_rows_spec = row_spec,
    n_col = n_cols_unit[[unit_idx]],
    case_col = case_cols_unit[[unit_idx]],
    ctrl_col = ctrl_cols_unit[[unit_idx]],
    effect_col = if (has_or_mode) or_cols_unit[[unit_idx]] else beta_cols_unit[[unit_idx]],
    se_col = if (has_or_mode) NA_integer_ else se_cols_unit[[unit_idx]],
    p_col = p_cols_unit[[unit_idx]],
    case = case_value,
    ctrl = ctrl_value,
    N = n_value,
    OR = or_value,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    OR_95CI = if (has_or_mode) effect_text else sprintf("%.2f(%.2f-%.2f)", or_value, ci_lower, ci_upper),
    beta = beta,
    se = se,
    p = p_value,
    source_p = p_raw,
    check.names = FALSE
  )
}

row_studies <- do.call(rbind, study_dat)
if (any(!is.finite(row_studies$beta)) || any(!is.finite(row_studies$se)) || any(row_studies$se <= 0)) {
  stop("selected rows ", row_spec, " have invalid beta/se values", call. = FALSE)
}
weights <- 1 / (row_studies$se ^ 2)
beta_meta <- sum(weights * row_studies$beta) / sum(weights)
se_meta <- sqrt(1 / sum(weights))
z_meta <- beta_meta / se_meta
p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)
ci_beta_low <- beta_meta - qnorm(0.975) * se_meta
ci_beta_high <- beta_meta + qnorm(0.975) * se_meta
q <- sum(weights * (row_studies$beta - beta_meta) ^ 2)
q_p <- pchisq(q, df = n_units - 1, lower.tail = FALSE)
i2 <- if (q <= (n_units - 1)) 0 else (q - (n_units - 1)) / q * 100

row_studies$record_type <- "input"
row_studies$meta_beta <- NA_real_
row_studies$meta_se <- NA_real_
row_studies$meta_z <- NA_real_
row_studies$meta_p <- NA_real_
row_studies$meta_ci_lower_beta <- NA_real_
row_studies$meta_ci_upper_beta <- NA_real_
row_studies$meta_OR <- NA_real_
row_studies$meta_OR_95CI <- NA_character_
row_studies$Q <- NA_real_
row_studies$Q_p <- NA_real_
row_studies$I2_percent <- NA_real_

meta_row <- row_studies[1, , drop = FALSE]
meta_row$record_type <- "meta"
meta_row$file_index <- NA_integer_
meta_row$study <- "meta-analysis"
meta_row$file <- paste(files, collapse = ";")
meta_row$source_row <- NA_integer_
meta_row$source_rows_spec <- row_spec
meta_row$n_col <- NA_integer_
meta_row$case_col <- NA_integer_
meta_row$ctrl_col <- NA_integer_
meta_row$effect_col <- NA_integer_
meta_row$se_col <- NA_integer_
meta_row$p_col <- NA_integer_
meta_row$case <- if (has_case_ctrl_mode) sum(row_studies$case, na.rm = TRUE) else NA_real_
meta_row$ctrl <- if (has_case_ctrl_mode) sum(row_studies$ctrl, na.rm = TRUE) else NA_real_
meta_row$N <- sum(row_studies$N, na.rm = TRUE)
meta_row$OR <- exp(beta_meta)
meta_row$CI_lower <- exp(ci_beta_low)
meta_row$CI_upper <- exp(ci_beta_high)
meta_row$OR_95CI <- sprintf("%.2f(%.2f-%.2f)", exp(beta_meta), exp(ci_beta_low), exp(ci_beta_high))
meta_row$beta <- beta_meta
meta_row$se <- se_meta
meta_row$p <- p_meta
meta_row$source_p <- format_meta_p_value(p_meta)
meta_row$meta_beta <- beta_meta
meta_row$meta_se <- se_meta
meta_row$meta_z <- z_meta
meta_row$meta_p <- p_meta
meta_row$meta_ci_lower_beta <- ci_beta_low
meta_row$meta_ci_upper_beta <- ci_beta_high
meta_row$meta_OR <- exp(beta_meta)
meta_row$meta_OR_95CI <- meta_row$OR_95CI
meta_row$Q <- q
meta_row$Q_p <- q_p
meta_row$I2_percent <- i2

result <- rbind(row_studies, meta_row)
result$record_type <- factor(result$record_type, levels = c("input", "meta"))
result <- result[order(result$meta_row, result$record_type, result$file_index), ]
result$record_type <- as.character(result$record_type)

keep_values <- list()
file1 <- tables[[1]]
file1_row <- unit_rows[[1]]
for (col in keep_cols) {
  values <- extract_col(file1, col, files[[1]])
  name <- names(file1)[[col]]
  keep_values[[name]] <- as.character(values[[file1_row]])
}
keep_output <- as.data.frame(keep_values, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(keep_output) > 0) {
  keep_output <- keep_output[rep(1, nrow(result)), , drop = FALSE]
  rownames(keep_output) <- NULL
}

if (full_output) {
  numeric_cols <- vapply(result, is.numeric, logical(1))
  result[numeric_cols] <- lapply(result[numeric_cols], function(x) ifelse(is.na(x), NA, signif(x, 8)))
  output <- cbind(result, keep_output)
} else {
  sample_output <- if (has_case_ctrl_mode) {
    data.frame(
      case = sprintf("%.0f", result$case),
      ctrl = sprintf("%.0f", result$ctrl),
      check.names = FALSE
    )
  } else {
    data.frame(N = sprintf("%.0f", result$N), check.names = FALSE)
  }
  output <- cbind(
    data.frame(study = result$study, check.names = FALSE),
    keep_output,
    sample_output,
    data.frame(
      `OR (95% CI)` = result$OR_95CI,
      `P-value` = ifelse(result$record_type == "input", result$source_p, format_meta_p_value(result$p)),
      check.names = FALSE
    )
  )
}

write.table(output, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote: ", normalizePath(out_file, mustWork = FALSE))

