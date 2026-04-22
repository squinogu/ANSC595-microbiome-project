
# =========================
# 02_random_forest.R
# =========================

library(phyloseq)
library(tidyverse)
library(caret)
library(randomForest)
library(pheatmap)
library(tibble)

# -------------------------
# 1. Set project folder
# -------------------------
project_dir <- "C:/Users/squinogu/my_project_qiime"
setwd(project_dir)

# -------------------------
# 2. Create output folders
# -------------------------
dir.create("results/random_forest", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/random_forest", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 3. Load phyloseq object
# -------------------------
ps <- readRDS("rds/ps_raw.rds")

meta <- data.frame(sample_data(ps))
if (!"treatment" %in% colnames(meta)) {
  stop("The metadata column 'treatment' was not found. Check your metadata.")
}

# -------------------------
# 4. Aggregate to genus level
# -------------------------
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 10, ps_genus)

# -------------------------
# 5. Create readable taxon names
# -------------------------
tax_df <- as.data.frame(tax_table(ps_genus))
tax_df[is.na(tax_df)] <- ""

make_tax_labels <- function(tax_df) {
  labels <- ifelse(
    tax_df$Genus != "",
    tax_df$Genus,
    ifelse(
      tax_df$Family != "",
      paste0("Unclassified_", tax_df$Family),
      ifelse(
        tax_df$Order != "",
        paste0("Unclassified_", tax_df$Order),
        paste0("Taxon_", seq_len(nrow(tax_df)))
      )
    )
  )
  make.unique(labels)
}

tax_labels <- make_tax_labels(tax_df)
taxa_names(ps_genus) <- tax_labels
rownames(tax_df) <- tax_labels

# -------------------------
# 6. Extract abundance table
# -------------------------
otu <- if (taxa_are_rows(ps_genus)) {
  t(as(otu_table(ps_genus), "matrix"))
} else {
  as(otu_table(ps_genus), "matrix")
}

otu <- otu[rowSums(otu) > 0, , drop = FALSE]

# -------------------------
# 7. Relative abundance + CLR transform
# -------------------------
otu_rel <- sweep(otu, 1, rowSums(otu), "/")
otu_rel[is.na(otu_rel)] <- 0

otu_clr <- log(otu_rel + 1e-6)
otu_clr <- otu_clr - rowMeans(otu_clr)

# -------------------------
# 8. Build modeling table
# -------------------------
meta <- data.frame(sample_data(ps_genus))
meta <- meta[rownames(otu_clr), , drop = FALSE]

rf_df <- as.data.frame(otu_clr, check.names = FALSE)
rf_df$treatment <- factor(meta$treatment)

# Put treatment first
rf_df <- rf_df[, c("treatment", setdiff(colnames(rf_df), "treatment"))]

# -------------------------
# 9. Cross-validated random forest
# -------------------------
min_n <- min(table(rf_df$treatment))
k_folds <- min(5, min_n)

if (k_folds < 2) {
  stop("You need at least 2 samples per treatment for cross-validation.")
}

set.seed(123)

ctrl <- trainControl(
  method = "cv",
  number = k_folds,
  savePredictions = "final"
)

rf_fit <- caret::train(
  treatment ~ .,
  data = rf_df,
  method = "rf",
  trControl = ctrl,
  importance = TRUE
)

saveRDS(rf_fit, "results/random_forest/rf_model.rds")

# -------------------------
# 10. Confusion matrix
# -------------------------
pred_df <- rf_fit$pred

if ("mtry" %in% colnames(pred_df)) {
  pred_df <- pred_df[pred_df$mtry == rf_fit$bestTune$mtry, , drop = FALSE]
}

cm <- confusionMatrix(pred_df$pred, pred_df$obs)

write.csv(as.data.frame(cm$table),
          "results/random_forest/confusion_matrix.csv",
          row.names = FALSE)

capture.output(cm,
               file = "results/random_forest/confusion_matrix_summary.txt")

write.csv(pred_df,
          "results/random_forest/cv_predictions.csv",
          row.names = FALSE)

acc_df <- data.frame(
  Accuracy = unname(cm$overall["Accuracy"]),
  Kappa = unname(cm$overall["Kappa"]),
  Folds = k_folds
)

write.csv(acc_df,
          "results/random_forest/model_accuracy.csv",
          row.names = FALSE)

# -------------------------
# 11. Variable importance
# -------------------------
imp <- as.data.frame(varImp(rf_fit)$importance)

if (ncol(imp) == 0) {
  stop("No variable importance values were returned.")
}

# Multiclass models often do not return 'Overall'
if (!"Overall" %in% colnames(imp)) {
  imp$Overall <- rowMeans(imp, na.rm = TRUE)
}

imp <- imp %>%
  tibble::rownames_to_column("Taxon") %>%
  arrange(desc(Overall))

write.csv(imp,
          "results/random_forest/variable_importance.csv",
          row.names = FALSE)

# -------------------------
# 12. Top 30 importance barplot
# -------------------------
top30_imp <- head(imp, 30)

p_bar <- ggplot(top30_imp, aes(x = reorder(Taxon, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Random Forest: Top 30 Important Taxa",
    x = NULL,
    y = "Importance"
  )

ggsave(
  filename = "figures/random_forest/top30_importance_barplot.png",
  plot = p_bar,
  width = 10,
  height = 8,
  dpi = 300
)
# -------------------------
# 13. Top 30 heatmap
# -------------------------
top30_taxa <- head(imp$Taxon, 30)

# Keep only taxa that really exist in rf_df
top30_taxa <- intersect(top30_taxa, colnames(rf_df))

if (length(top30_taxa) == 0) {
  stop("None of the top taxa were found in rf_df column names.")
}

heat_mat <- t(as.matrix(rf_df[, top30_taxa, drop = FALSE]))

ann_col <- data.frame(Treatment = rf_df$treatment)
rownames(ann_col) <- rownames(rf_df)

pheatmap(
  mat = heat_mat,
  annotation_col = ann_col,
  scale = "row",
  show_colnames = FALSE,
  fontsize_row = 8,
  filename = "figures/random_forest/top30_heatmap.png",
  width = 10,
  height = 8
)

cat("Done. Random Forest results saved in results/random_forest and figures/random_forest\n")