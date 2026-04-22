install.packages(c(
  "tidyverse", "randomForest", "caret", "pheatmap",
  "igraph", "Hmisc", "remotes"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("phyloseq", "DESeq2"), ask = FALSE)

remotes::install_github("jbisanz/qiime2R")

library(qiime2R)
library(phyloseq)
library(tidyverse)

# =========================
project_dir   <- "C:/Users/squinogu/my_project_qiime"
metadata_file <- "mymetadata.tsv"
table_file    <- "table.qza"
taxonomy_file <- "taxonomy.qza"
tree_file     <- "rooted-tree.qza"
# =========================

setwd(project_dir)

dir.create("scripts", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("rds", showWarnings = FALSE)

ps <- qza_to_phyloseq(
  features = table_file,
  taxonomy = taxonomy_file,
  tree = tree_file,
  metadata = metadata_file
)

meta <- data.frame(sample_data(ps), check.names = FALSE)

# Change this if your metadata column is NOT called treatment
stopifnot("treatment" %in% colnames(meta))

meta$treatment <- as.character(meta$treatment)
meta$treatment <- gsub("\\+", "_plus_", meta$treatment)

# Edit levels only if your treatment names are different
meta$treatment <- factor(
  meta$treatment,
  levels = c("Blank", "BP", "FMA", "FMA_plus_BP", "FMB", "FMB_plus_BP")
)

sample_data(ps) <- sample_data(meta)

tax <- as.data.frame(tax_table(ps))
tax[is.na(tax)] <- ""
for (i in seq_len(ncol(tax))) {
  tax[, i] <- gsub("^[A-Za-z]__", "", tax[, i])
}
tax_table(ps) <- as.matrix(tax)

saveRDS(ps, "rds/ps_raw.rds")
write.csv(meta, "results/checked_metadata.csv", row.names = TRUE)

cat("Done. Saved phyloseq object to rds/ps_raw.rds\n")