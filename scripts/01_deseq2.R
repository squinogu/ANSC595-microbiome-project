
library(phyloseq)
library(DESeq2)
library(tidyverse)

setwd("C:/Users/squinogu/my_project_qiime")

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/deseq2", recursive = TRUE, showWarnings = FALSE)

ps <- readRDS("rds/ps_raw.rds")

make_tax_labels <- function(tax_df) {
  out <- ifelse(
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
  make.unique(out)
}

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 10, ps_genus)

tax_df <- as.data.frame(tax_table(ps_genus))
tax_labels <- make_tax_labels(tax_df)
taxa_names(ps_genus) <- tax_labels
rownames(tax_df) <- tax_labels

dds <- phyloseq_to_deseq2(ps_genus, ~ treatment)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

treat_levels <- levels(colData(dds)$treatment)
pair_list <- combn(treat_levels, 2, simplify = FALSE)

summary_list <- list()

for (pair in pair_list) {
  g1 <- pair[1]
  g2 <- pair[2]
  
  res <- results(
    dds,
    contrast = c("treatment", g1, g2),
    alpha = 0.05,
    cooksCutoff = FALSE
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Taxon") %>%
    left_join(
      tax_df %>% rownames_to_column("Taxon"),
      by = "Taxon"
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  out_name <- paste0(g1, "_vs_", g2)
  
  write.csv(
    res_df,
    file = file.path("results/deseq2", paste0(out_name, ".csv")),
    row.names = FALSE
  )
  
  sig_df <- res_df %>%
    filter(!is.na(padj) & padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange)))
  
  summary_list[[out_name]] <- data.frame(
    comparison = out_name,
    n_sig_taxa = nrow(sig_df)
  )
  
  if (nrow(sig_df) > 0) {
    plot_df <- sig_df %>%
      slice_head(n = 15) %>%
      mutate(
        Label = ifelse(Genus != "", Genus, Taxon),
        Direction = ifelse(log2FoldChange > 0,
                           paste0("Enriched in ", g1),
                           paste0("Enriched in ", g2))
      )
    
    p <- ggplot(plot_df, aes(x = reorder(Label, log2FoldChange),
                             y = log2FoldChange,
                             fill = Direction)) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(
        title = paste(g1, "vs", g2),
        x = NULL,
        y = "log2 fold change"
      )
    
    ggsave(
      filename = file.path("figures/deseq2", paste0(out_name, "_top15.png")),
      plot = p,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}

bind_rows(summary_list) %>%
  write.csv("results/deseq2/deseq2_summary.csv", row.names = FALSE)

cat("Done. DESeq2 results saved in results/deseq2 and figures/deseq2\n")
