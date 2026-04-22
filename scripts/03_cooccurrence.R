# =========================
# 03_cooccurrence.R
# =========================

library(phyloseq)
library(tidyverse)
library(igraph)
library(tibble)

# -------------------------
# 1. Set project folder
# -------------------------
project_dir <- "C:/Users/squinogu/my_project_qiime"
setwd(project_dir)

# -------------------------
# 2. Create output folders
# -------------------------
dir.create("results/cooccurrence", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/cooccurrence", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 3. Load phyloseq object
# -------------------------
ps <- readRDS("rds/ps_raw.rds")

# -------------------------
# 4. Aggregate to genus level
# -------------------------
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 20, ps_genus)

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
# 7. Convert to relative abundance
# -------------------------
otu_rel <- sweep(otu, 1, rowSums(otu), "/")
otu_rel[is.na(otu_rel)] <- 0

# -------------------------
# 8. Filter taxa for network
#    prevalence >= 20%
#    mean relative abundance >= 0.001
# -------------------------
prev <- colMeans(otu_rel > 0)
mean_abund <- colMeans(otu_rel)

keep_taxa <- names(prev[prev >= 0.20 & mean_abund >= 0.001])

otu_f <- otu_rel[, keep_taxa, drop = FALSE]

if (ncol(otu_f) < 3) {
  stop("Too few taxa remained after filtering. Lower the filtering thresholds.")
}

# -------------------------
# 9. Spearman correlation matrix
# -------------------------
rho_mat <- cor(otu_f, method = "spearman")

# -------------------------
# 10. P-value matrix with pairwise cor.test
# -------------------------
n_taxa <- ncol(otu_f)

p_mat <- matrix(NA, nrow = n_taxa, ncol = n_taxa,
                dimnames = list(colnames(otu_f), colnames(otu_f)))

for (i in 1:(n_taxa - 1)) {
  for (j in (i + 1):n_taxa) {
    test_out <- suppressWarnings(cor.test(
      otu_f[, i],
      otu_f[, j],
      method = "spearman",
      exact = FALSE
    ))
    p_mat[i, j] <- test_out$p.value
    p_mat[j, i] <- test_out$p.value
  }
}

diag(p_mat) <- 0

# BH adjustment
upper_idx <- upper.tri(p_mat)
padj_vals <- p.adjust(p_mat[upper_idx], method = "BH")

padj_mat <- matrix(1, nrow = n_taxa, ncol = n_taxa,
                   dimnames = list(colnames(otu_f), colnames(otu_f)))
padj_mat[upper_idx] <- padj_vals
padj_mat[lower.tri(padj_mat)] <- t(padj_mat)[lower.tri(padj_mat)]
diag(padj_mat) <- 0

# -------------------------
# 11. Build edge list
#    first try |rho| >= 0.6
#    if empty, try |rho| >= 0.5
# -------------------------
build_edges <- function(rho_threshold) {
  idx <- which(abs(rho_mat) >= rho_threshold & padj_mat < 0.05, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  
  if (nrow(idx) == 0) {
    return(data.frame())
  }
  
  data.frame(
    from = rownames(rho_mat)[idx[, 1]],
    to = colnames(rho_mat)[idx[, 2]],
    rho = rho_mat[idx],
    padj = padj_mat[idx],
    sign = ifelse(rho_mat[idx] > 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )
}

edges <- build_edges(0.6)
rho_used <- 0.6

if (nrow(edges) == 0) {
  edges <- build_edges(0.5)
  rho_used <- 0.5
}

if (nrow(edges) == 0) {
  stop("No significant edges were found even at |rho| >= 0.5. We may need to relax the thresholds.")
}

# -------------------------
# 12. Build node table
# -------------------------
nodes <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE) %>%
  left_join(tax_df %>% tibble::rownames_to_column("name"), by = "name")

g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

nodes$degree <- degree(g)[nodes$name]
nodes$betweenness <- betweenness(g)[nodes$name]
nodes$closeness <- closeness(g, normalized = TRUE)[nodes$name]

# -------------------------
# 13. Save tables
# -------------------------
write.csv(edges, "results/cooccurrence/cooccurrence_edges.csv", row.names = FALSE)
write.csv(nodes, "results/cooccurrence/cooccurrence_nodes.csv", row.names = FALSE)

summary_lines <- c(
  paste("Number of taxa in filtered table:", ncol(otu_f)),
  paste("Number of samples:", nrow(otu_f)),
  paste("Correlation threshold used:", rho_used),
  paste("Number of edges:", nrow(edges)),
  paste("Number of nodes:", nrow(nodes)),
  paste("Positive edges:", sum(edges$sign == 'positive')),
  paste("Negative edges:", sum(edges$sign == 'negative'))
)

writeLines(summary_lines, "results/cooccurrence/cooccurrence_summary.txt")

# -------------------------
# 14. Network plot
# -------------------------
set.seed(123)
lay <- layout_with_fr(g)

png("figures/cooccurrence/cooccurrence_network.png",
    width = 2000, height = 1600, res = 220)

plot(
  g,
  layout = lay,
  vertex.size = 5 + nodes$degree[match(V(g)$name, nodes$name)] * 1.5,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  vertex.frame.color = "gray40",
  edge.width = abs(E(g)$rho) * 4,
  edge.color = ifelse(E(g)$sign == "positive", "steelblue", "tomato"),
  main = paste0("Co-occurrence network (|rho| >= ", rho_used, ", BH-adjusted p < 0.05)")
)

dev.off()

# -------------------------
# 15. Top 20 nodes by degree
# -------------------------
top_nodes <- nodes %>%
  arrange(desc(degree)) %>%
  slice_head(n = 20)

p_top <- ggplot(top_nodes, aes(x = reorder(name, degree), y = degree)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top 20 taxa by network degree",
    x = NULL,
    y = "Degree"
  )

ggsave(
  filename = "figures/cooccurrence/top20_nodes_degree.png",
  plot = p_top,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Done. Co-occurrence results saved in results/cooccurrence and figures/cooccurrence\n")