# ============================================================
# Prepare MaxQuant proteinGroups.txt for NAguideR
# MCM-1G versus MCM-DLN
# ============================================================

# -----------------------------
# 1. Working directory and files
# -----------------------------

setwd("~/Proteomics_Melanoma_Pauser")

input_file <- "proteinGroups.txt"
output_dir <- "NAguideR_input"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------
# 2. Read MaxQuant proteinGroups.txt
# -----------------------------

protein_groups <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

cat(
  "Raw proteinGroups dimensions:",
  nrow(protein_groups), "proteins x",
  ncol(protein_groups), "columns\n"
)

# -----------------------------
# 3. Remove reverse hits and contaminants
# -----------------------------

required_annotation_columns <- c(
  "Majority protein IDs",
  "Reverse",
  "Potential contaminant"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(protein_groups)
)

if (length(missing_annotation_columns) > 0) {
  stop(
    "The following required MaxQuant columns are missing: ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

protein_groups <- protein_groups[
  (is.na(protein_groups$Reverse) |
     protein_groups$Reverse != "+") &
    (is.na(protein_groups$`Potential contaminant`) |
       protein_groups$`Potential contaminant` != "+"),
  ,
  drop = FALSE
]

cat(
  "Proteins after removal of contaminants and reverse hits:",
  nrow(protein_groups), "\n"
)

# -----------------------------
# 4. Extract melanoma LFQ columns
# -----------------------------

# Expected sample names:
# 1G_A_01, 1G_A_02, 1G_A_03,
# 1G_B_01, ..., 1G_C_03,
# DLN_A_01, ..., DLN_C_03

lfq_pattern <- paste0(
  "^LFQ intensity ",
  "(1G_A|1G_B|1G_C|DLN_A|DLN_B|DLN_C)_",
  "(01|02|03)$"
)

lfq_columns <- grep(
  lfq_pattern,
  colnames(protein_groups),
  value = TRUE
)

if (length(lfq_columns) != 18) {
  stop(
    "Expected 18 LFQ intensity columns, but found ",
    length(lfq_columns),
    ".\nDetected columns:\n",
    paste(lfq_columns, collapse = "\n")
  )
}

lfq_raw <- protein_groups[, lfq_columns, drop = FALSE]

# Use protein IDs as row names
protein_ids <- protein_groups$`Majority protein IDs`

if (any(is.na(protein_ids) | protein_ids == "")) {
  stop("Missing entries were detected in the Majority protein IDs column.")
}

rownames(lfq_raw) <- make.unique(protein_ids)

# -----------------------------
# 5. Clean sample names
# -----------------------------

colnames(lfq_raw) <- sub(
  "^LFQ intensity ",
  "",
  colnames(lfq_raw)
)

# Arrange samples in the intended order
desired_order <- c(
  "1G_A_01", "1G_A_02", "1G_A_03",
  "1G_B_01", "1G_B_02", "1G_B_03",
  "1G_C_01", "1G_C_02", "1G_C_03",
  "DLN_A_01", "DLN_A_02", "DLN_A_03",
  "DLN_B_01", "DLN_B_02", "DLN_B_03",
  "DLN_C_01", "DLN_C_02", "DLN_C_03"
)

missing_samples <- setdiff(
  desired_order,
  colnames(lfq_raw)
)

if (length(missing_samples) > 0) {
  stop(
    "The following expected samples were not found: ",
    paste(missing_samples, collapse = ", ")
  )
}

lfq_raw <- lfq_raw[, desired_order, drop = FALSE]

# -----------------------------
# 6. Convert LFQ values to numeric
# -----------------------------

lfq <- as.data.frame(
  lapply(
    lfq_raw,
    function(x) suppressWarnings(as.numeric(x))
  ),
  check.names = FALSE
)

rownames(lfq) <- rownames(lfq_raw)

# MaxQuant uses zero to indicate a missing LFQ value
lfq[lfq == 0] <- NA

# -----------------------------
# 7. Filter proteins
# -----------------------------

# Require at least two valid LFQ measurements in each condition.
# Each condition contains nine measurement columns.

cols_1G <- grep("^1G_", colnames(lfq))
cols_DLN <- grep("^DLN_", colnames(lfq))

valid_1G <- rowSums(!is.na(lfq[, cols_1G, drop = FALSE]))
valid_DLN <- rowSums(!is.na(lfq[, cols_DLN, drop = FALSE]))

keep_proteins <- valid_1G >= 2 & valid_DLN >= 2

lfq_filtered <- lfq[keep_proteins, , drop = FALSE]

cat(
  "Proteins retained after filtering:",
  nrow(lfq_filtered), "\n"
)

# -----------------------------
# 8. Log2 transformation
# -----------------------------

# MaxQuant LFQ intensities are normally provided on a linear scale.
# Values are therefore log2-transformed before NAguideR analysis.

if (all(is.na(lfq_filtered))) {
  stop("The filtered LFQ matrix contains no valid values.")
}

if (max(as.matrix(lfq_filtered), na.rm = TRUE) > 1000) {
  lfq_log2 <- log2(lfq_filtered)
  transformation_status <- "LFQ intensities were log2-transformed."
} else {
  lfq_log2 <- lfq_filtered
  transformation_status <- paste(
    "The values did not appear to be on the raw LFQ intensity scale;",
    "no additional log2 transformation was applied."
  )
}

cat(transformation_status, "\n")

# -----------------------------
# 9. Remove proteins missing in all samples
# -----------------------------

lfq_log2 <- lfq_log2[
  rowSums(!is.na(lfq_log2)) > 0,
  ,
  drop = FALSE
]

# -----------------------------
# 10. Generate metadata
# -----------------------------

sample_names <- colnames(lfq_log2)

metadata_full <- data.frame(
  Sample = sample_names,
  Condition = ifelse(
    grepl("^1G_", sample_names),
    "MCM-1G",
    "MCM-DLN"
  ),
  Technical_replicate = sub(
    "^(1G|DLN)_([ABC])_[0-9]+$",
    "\\2",
    sample_names
  ),
  Measurement_replicate = sub(
    "^.*_([0-9]+)$",
    "\\1",
    sample_names
  ),
  stringsAsFactors = FALSE
)

# Metadata format required for NAguideR
metadata_NAguideR <- data.frame(
  Samples = metadata_full$Sample,
  Groups = metadata_full$Condition,
  stringsAsFactors = FALSE
)

# Confirm exact agreement between input matrix and metadata
stopifnot(
  identical(
    colnames(lfq_log2),
    metadata_NAguideR$Samples
  )
)

# -----------------------------
# 11. Final quality checks
# -----------------------------

stopifnot(is.data.frame(lfq_log2))
stopifnot(all(vapply(lfq_log2, is.numeric, logical(1))))
stopifnot(ncol(lfq_log2) == 18)
stopifnot(sum(lfq_log2 == 0, na.rm = TRUE) == 0)
stopifnot(!anyDuplicated(rownames(lfq_log2)))
stopifnot(!anyDuplicated(metadata_NAguideR$Samples))

missing_values_per_sample <- colSums(is.na(lfq_log2))
missing_percent_per_sample <- round(
  colMeans(is.na(lfq_log2)) * 100,
  2
)

missing_summary <- data.frame(
  Sample = colnames(lfq_log2),
  Missing_values = missing_values_per_sample,
  Missing_percent = missing_percent_per_sample,
  stringsAsFactors = FALSE
)

overall_summary <- data.frame(
  Proteins = nrow(lfq_log2),
  Samples = ncol(lfq_log2),
  Total_missing_values = sum(is.na(lfq_log2)),
  Overall_missing_percent = round(
    mean(is.na(lfq_log2)) * 100,
    2
  )
)

# -----------------------------
# 12. Save NAguideR input files
# -----------------------------

write.csv(
  lfq_log2,
  file = file.path(
    output_dir,
    "NAguideR_input_Melanoma_LFQ_log2_NA.csv"
  ),
  row.names = TRUE,
  na = ""
)

write.csv(
  metadata_NAguideR,
  file = file.path(
    output_dir,
    "NAguideR_metadata_Melanoma_STRICT.csv"
  ),
  row.names = FALSE
)

write.csv(
  metadata_full,
  file = file.path(
    output_dir,
    "Melanoma_sample_metadata_full.csv"
  ),
  row.names = FALSE
)

write.csv(
  missing_summary,
  file = file.path(
    output_dir,
    "Missing_values_per_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  overall_summary,
  file = file.path(
    output_dir,
    "Missing_values_overall.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# 13. Print summary
# -----------------------------

cat("\nNAguideR files successfully generated.\n\n")

cat(
  "LFQ input file:\n",
  file.path(
    output_dir,
    "NAguideR_input_Melanoma_LFQ_log2_NA.csv"
  ),
  "\n\n"
)

cat(
  "Metadata file:\n",
  file.path(
    output_dir,
    "NAguideR_metadata_Melanoma_STRICT.csv"
  ),
  "\n\n"
)

cat(
  "Final dimensions:",
  nrow(lfq_log2), "proteins x",
  ncol(lfq_log2), "measurement replicates\n"
)

cat(
  "Total missing values:",
  sum(is.na(lfq_log2)), "\n"
)

cat(
  "Overall missingness:",
  round(mean(is.na(lfq_log2)) * 100, 2),
  "%\n"
)

print(metadata_full)
print(missing_summary)
print(overall_summary)


# Launch the graphical interface
NAguideR::NAguideR_app()