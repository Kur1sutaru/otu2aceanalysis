# ============================================================
# ANÁLISE DE MICROBIOTA INTESTINAL — CCL vs CONTROLE
# Sequenciamento 16S rRNA (IonTorrent) | Região V4
# ============================================================
# Autora: [seu nome]
# Data:   2026
# R:      4.3.2
#
# Pacotes necessários:
#   install.packages(c("phyloseq","vegan","pheatmap","readxl",
#                      "dplyr","tibble","ggplot2","RColorBrewer"))
#   BiocManager::install(c("phyloseq","DESeq2","microbiomeMarker"))
#
# Arquivos de entrada (D:/metagenomics/):
#   - ASVs_counts_ccl.tsv
#   - ASVs_counts_control.tsv
#   - ASVs_taxonomy_CCL.csv
#   - ASVs_taxonomy_Controles.csv
#   - DADOS SOCIODEMOGRÁFICOS - CCL.xlsx
#   - DADOS SOCIODEMOGRÁFICOS - CONTROLE.xlsx
# ============================================================

library(phyloseq)
library(vegan)
library(pheatmap)
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

# ── Configurações ─────────────────────────────────────────────
setwd("D:/metagenomics")
set.seed(42)
TOP_N  <- 30    # OTUs para heatmap
OUTDIR <- "D:/metagenomics/outputs"
dir.create(OUTDIR, showWarnings = FALSE)

# ── Grupos corretos (fonte da verdade = tabelas clínicas) ─────
# ATENÇÃO: o pipeline separou incorretamente algumas amostras.
# A correção é feita aqui, sem alterar os dados de sequenciamento.

CCL_IDS <- sprintf("IonXpress_%03d", c(3,4,5,8,9,10,22,23,24))
CTRL_IDS <- sprintf("IonXpress_%03d",
             c(1,2,6,7,11,12,13,14,15,16,17,18,19,20,21,25,26))

# ============================================================
# PASSO 1: CARREGAR DADOS E CONSTRUIR OBJETO PHYLOSEQ
# ============================================================
message("== PASSO 1: Carregando dados ==")

# Contagens
counts_ccl  <- read.delim("ASVs_counts_ccl.tsv",     row.names = 1, check.names = FALSE)
counts_ctrl <- read.delim("ASVs_counts_control.tsv",  row.names = 1, check.names = FALSE)

# Taxonomia (posição = IDs de ASV)
tax_ccl  <- read.csv("ASVs_taxonomy_CCL.csv",       stringsAsFactors = FALSE)
tax_ctrl <- read.csv("ASVs_taxonomy_Controles.csv",  stringsAsFactors = FALSE)
rownames(tax_ccl)  <- rownames(counts_ccl)
rownames(tax_ctrl) <- rownames(counts_ctrl)

# Renomear colunas para IonXpress_XXX
rename_to_ionxpress <- function(col) {
  num <- as.integer(gsub("\\D", "", strsplit(col, "_|-")[[1]][1]))
  sprintf("IonXpress_%03d", num)
}
colnames(counts_ccl)  <- sapply(colnames(counts_ccl),  rename_to_ionxpress)
colnames(counts_ctrl) <- sapply(colnames(counts_ctrl), rename_to_ionxpress)

# Combinar matrizes (ASVs não compartilhados recebem 0)
all_asvs <- union(rownames(counts_ccl), rownames(counts_ctrl))

mat_ccl  <- matrix(0L, nrow = length(all_asvs), ncol = ncol(counts_ccl),
                    dimnames = list(all_asvs, colnames(counts_ccl)))
mat_ctrl <- matrix(0L, nrow = length(all_asvs), ncol = ncol(counts_ctrl),
                    dimnames = list(all_asvs, colnames(counts_ctrl)))

mat_ccl[rownames(counts_ccl),  ] <- as.matrix(counts_ccl)
mat_ctrl[rownames(counts_ctrl), ] <- as.matrix(counts_ctrl)

counts_all <- cbind(mat_ccl, mat_ctrl)
counts_all <- counts_all[, colnames(counts_all) != "IonXpress_027"]  # sem dados clínicos

# Taxonomia combinada
tax_combined <- rbind(
  tax_ccl[rownames(tax_ccl) %in% rownames(counts_ccl), ],
  tax_ctrl[!rownames(tax_ctrl) %in% rownames(counts_ccl), ]
)
tax_combined <- tax_combined[rownames(counts_all), , drop = FALSE]

# Grupo correto por amostra
grupo_vec <- ifelse(colnames(counts_all) %in% CCL_IDS, "CCL", "Controle")
sdata <- data.frame(Grupo = factor(grupo_vec, levels = c("CCL","Controle")),
                    row.names = colnames(counts_all))

# Objeto phyloseq
ps <- phyloseq(
  otu_table(counts_all, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_combined)),
  sample_data(sdata)
)

cat(sprintf("Amostras: %d | ASVs: %d\nCCL: %d | Controle: %d\n",
            nsamples(ps), ntaxa(ps),
            sum(sdata$Grupo == "CCL"),
            sum(sdata$Grupo == "Controle")))

# ============================================================
# PASSO 2: RAREFAÇÃO E NORMALIZAÇÃO
# ============================================================
message("\n== PASSO 2: Rarefação ==")

min_depth <- min(sample_sums(ps))
cat(sprintf("Profundidade mínima: %d reads\n", min_depth))

ps_r   <- rarefy_even_depth(ps,
            sample.size = min_depth,
            replace     = FALSE,
            trimOTUs    = TRUE,
            verbose     = TRUE)

ps_rel <- transform_sample_counts(ps_r, function(x) x / sum(x))

# ============================================================
# PASSO 3: DIVERSIDADE ALFA
# ============================================================
message("\n== PASSO 3: Diversidade Alfa ==")

alpha_div <- estimate_richness(ps_r,
               measures = c("Observed", "Shannon", "Simpson"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div$Grupo    <- sample_data(ps_r)$Grupo

# Teste Mann-Whitney U por índice
for (idx in c("Observed", "Shannon", "Simpson")) {
  x_ccl  <- alpha_div[[idx]][alpha_div$Grupo == "CCL"]
  x_ctrl <- alpha_div[[idx]][alpha_div$Grupo == "Controle"]
  res    <- wilcox.test(x_ccl, x_ctrl, exact = FALSE)
  cat(sprintf("%s: CCL=%.3f±%.3f vs Ctrl=%.3f±%.3f | p=%.4f\n",
              idx,
              mean(x_ccl), sd(x_ccl),
              mean(x_ctrl), sd(x_ctrl),
              res$p.value))
}

# Boxplot alfa diversidade
df_alpha_long <- tidyr::pivot_longer(alpha_div,
                   cols = c("Observed","Shannon","Simpson"),
                   names_to = "Index", values_to = "Value")

p_alpha <- ggplot(df_alpha_long, aes(x = Grupo, y = Value, fill = Grupo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~Index, scales = "free_y") +
  scale_fill_manual(values = c(CCL = "#E07070", Controle = "#70B8D4")) +
  labs(title = "Diversidade Alfa — CCL vs Controle",
       x = NULL, y = "Valor do índice") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 11))

ggsave(file.path(OUTDIR, "AlphaDiversity_Boxplot.pdf"),
       p_alpha, width = 10, height = 4, dpi = 300)
ggsave(file.path(OUTDIR, "AlphaDiversity_Boxplot.png"),
       p_alpha, width = 10, height = 4, dpi = 300)

# Salvar tabela
write.xlsx2_safe <- function(df, file) {
  tryCatch(writexl::write_xlsx(df, file), error = function(e) write.csv(df, sub(".xlsx",".csv",file)))
}
write.csv(alpha_div, file.path(OUTDIR, "Tabela_AlphaDiversity_R.csv"), row.names = FALSE)
message("  AlphaDiversity_Boxplot.pdf/png + Tabela_AlphaDiversity_R.csv")

# ============================================================
# PASSO 4: DIVERSIDADE BETA (PCoA Bray-Curtis + PERMANOVA)
# ============================================================
message("\n== PASSO 4: Diversidade Beta ==")

# PCoA
otu_mat_r <- as(otu_table(ps_r), "matrix")
if (!taxa_are_rows(ps_r)) otu_mat_r <- t(otu_mat_r)

bc_dist <- vegdist(t(otu_mat_r), method = "bray")
pcoa    <- cmdscale(bc_dist, k = 2, eig = TRUE)

# Variância explicada
eig_vals  <- pcoa$eig[pcoa$eig > 0]
var_expl  <- round(eig_vals / sum(eig_vals) * 100, 1)

df_pcoa <- data.frame(
  PC1    = pcoa$points[, 1],
  PC2    = pcoa$points[, 2],
  Grupo  = sample_data(ps_r)$Grupo,
  Sample = rownames(pcoa$points)
)

p_pcoa <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Grupo)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(level = 0.95, lwd = 0.8) +
  scale_color_manual(values = c(CCL = "#E07070", Controle = "#70B8D4")) +
  labs(title = "Diversidade Beta — PCoA (Bray-Curtis)",
       x = sprintf("[%.1f%%]", var_expl[1]),
       y = sprintf("[%.1f%%]", var_expl[2])) +
  theme_classic() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(OUTDIR, "BetaDiversity_PCoA.pdf"),
       p_pcoa, width = 7, height = 5, dpi = 300)
ggsave(file.path(OUTDIR, "BetaDiversity_PCoA.png"),
       p_pcoa, width = 7, height = 5, dpi = 300)

# PERMANOVA
permanova <- adonis2(bc_dist ~ Grupo, data = as.data.frame(sample_data(ps_r)),
                     permutations = 999)
cat("\nPERMANOVA (Bray-Curtis ~ Grupo):\n")
print(permanova)

# PERMDISP (homogeneidade de dispersão)
dispr <- betadisper(bc_dist, sample_data(ps_r)$Grupo)
cat("\nPERMDISP p-valor:", permutest(dispr)$tab["Groups","Pr(>F)"], "\n")

message("  BetaDiversity_PCoA.pdf/png")

# ============================================================
# PASSO 5: HEATMAP (pheatmap) — Top OTUs × Amostras
# ============================================================
message("\n== PASSO 5: Heatmap Top OTUs ==")

# Top N OTUs
top_taxa  <- names(sort(taxa_sums(ps_rel), decreasing = TRUE)[1:TOP_N])
ps_top    <- prune_taxa(top_taxa, ps_rel)

otu_matrix <- as(otu_table(ps_top), "matrix")
if (!taxa_are_rows(ps_top)) otu_matrix <- t(otu_matrix)

# Labels de gênero
tax_df     <- as.data.frame(tax_table(ps_top))
tax_labels <- ifelse(
  !is.na(tax_df$Genus) & tax_df$Genus != "",
  paste0(tax_df$Genus, " (", rownames(tax_df), ")"),
  ifelse(!is.na(tax_df$Family) & tax_df$Family != "",
         paste0("[", tax_df$Family, "] ", rownames(tax_df)),
         rownames(tax_df))
)
rownames(otu_matrix) <- tax_labels

# Ordenar: Controle | CCL (clusterizado dentro de cada grupo)
ctrl_cols <- colnames(otu_matrix)[colnames(otu_matrix) %in% CTRL_IDS]
ccl_cols  <- colnames(otu_matrix)[colnames(otu_matrix) %in% CCL_IDS]
col_order <- c(ctrl_cols, ccl_cols)
otu_matrix <- otu_matrix[, col_order]

# Anotação de coluna
annotation_col <- data.frame(
  Grupo = ifelse(col_order %in% CCL_IDS, "CCL", "Controle"),
  row.names = col_order
)
ann_colors <- list(Grupo = c(CCL = "#E07070", Controle = "#70B8D4"))

# Heatmap com escala Z-score por linha
otu_log <- log10(otu_matrix * 1000 + 1)

pdf(file.path(OUTDIR, "Heatmap_R_TopOTUs.pdf"), width = 16, height = 11)
pheatmap(otu_log,
         annotation_col           = annotation_col,
         annotation_colors        = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method        = "ward.D2",
         scale        = "row",
         color        = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         border_color = NA,
         fontsize_row = 8, fontsize_col = 9,
         angle_col    = 45, cellwidth = 20, cellheight = 13,
         gaps_col     = length(ctrl_cols),   # linha vertical entre grupos
         main         = sprintf("Top %d OTUs | CCL vs Controle", TOP_N))
dev.off()

# PNG também
png(file.path(OUTDIR, "Heatmap_R_TopOTUs.png"),
    width = 16, height = 11, units = "in", res = 300)
pheatmap(otu_log,
         annotation_col = annotation_col, annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2", scale = "row",
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         border_color = NA, fontsize_row = 8, fontsize_col = 9,
         angle_col = 45, cellwidth = 20, cellheight = 13,
         gaps_col = length(ctrl_cols),
         main = sprintf("Top %d OTUs | CCL vs Controle", TOP_N))
dev.off()

message("  Heatmap_R_TopOTUs.pdf/png")

# ============================================================
# PASSO 6: CORRELAÇÃO SPEARMAN — OTUs × ACE + Sociodemográfico
# ============================================================
message("\n== PASSO 6: Correlação Spearman ==")

# Dados ACE (subdomínios + total)
ace_cols_vec <- c("Orientação","Atenção","Recordação","Mem_anterógrada",
                  "Mem_retrógrada","Rec_reconhecimento","Fluência_verbal",
                  "Compreensão","Escrita","Repetição","Nomeação","Leitura",
                  "Visual_espacial","ACE_Total")

ace_ctrl_data <- list(
  IonXpress_001=c(13,5,2,7,4,9,12,6,1,4,12,1,8,92), IonXpress_002=c(13,5,0,7,4,11,13,7,1,4,12,1,8,94),
  IonXpress_015=c(10,4,3,7,4,8,7,6,1,4,11,1,1,72),  IonXpress_019=c(13,5,3,7,4,7,13,8,1,4,12,1,8,94),
  IonXpress_006=c(13,4,2,7,3,10,12,8,1,4,12,1,8,93), IonXpress_007=c(12,5,3,5,3,11,7,7,1,3,10,0,6,81),
  IonXpress_026=c(13,4,2,7,4,12,13,8,1,4,12,1,8,97), IonXpress_018=c(13,1,2,4,3,3,10,7,1,2,11,1,2,68),
  IonXpress_020=c(13,1,3,6,4,5,14,8,1,3,12,1,6,85),  IonXpress_013=c(13,5,1,6,3,8,6,4,1,4,12,1,8,80),
  IonXpress_016=c(13,4,3,6,4,6,12,8,1,3,12,1,8,89),  IonXpress_017=c(13,5,3,7,4,10,7,4,1,4,12,1,8,87),
  IonXpress_014=c(13,5,1,7,4,11,11,6,1,4,12,1,8,92), IonXpress_011=c(13,5,2,5,4,7,12,8,1,4,11,1,7,88),
  IonXpress_025=c(13,5,3,6,4,9,14,8,1,4,12,1,6,94),  IonXpress_012=c(13,5,2,7,4,9,9,8,1,4,12,1,8,91),
  IonXpress_021=c(13,3,3,7,1,9,12,8,1,3,10,1,7,86)
)
ace_ccl_data <- list(
  IonXpress_009=c(13,4,2,6,1,5,8,5,1,4,10,1,6,72),   IonXpress_010=c(13,3,2,7,3,6,7,7,1,4,12,1,8,81),
  IonXpress_023=c(12,2,1,5,3,5,4,3,1,1,7,0,7,59),    IonXpress_003=c(13,4,0,5,3,3,8,6,1,3,9,1,4,67),
  IonXpress_004=c(12,5,0,5,2,4,2,7,1,2,11,1,2,62),   IonXpress_005=c(12,0,2,5,1,5,2,7,1,2,10,1,3,59),
  IonXpress_008=c(13,1,1,6,3,4,7,8,1,4,12,1,8,77),   IonXpress_024=c(13,5,2,7,4,9,12,8,1,4,11,1,1,86),
  IonXpress_022=c(11,5,3,0,1,2,13,8,1,4,11,1,8,76)
)

ace_all_data <- c(ace_ctrl_data, ace_ccl_data)
ace_df <- as.data.frame(t(do.call(rbind,
           lapply(ace_all_data, function(x) setNames(x, ace_cols_vec)))))
colnames(ace_df) <- names(ace_all_data)

# Sociodemográfico
load_socio_r <- function(filepath, grupo) {
  df <- read_excel(filepath, sheet = "Página1")
  df <- df[!is.na(df[[1]]) & !grepl("LEGENDA", df[[1]]), ]
  df$SampleID <- sprintf("IonXpress_%03d",
                  as.integer(gsub("\\D", "", df[[1]])))
  df$Grupo <- grupo
  df[, c("SampleID", "Grupo",
         "Idade (anos)", "Anos de estudo", "Doenças crônicas?",
         "Realiza atividade física?", "Consome bebidas alcoolicas")] |>
    setNames(c("SampleID","Grupo","Idade","Anos_estudo",
               "Doencas","Ativ_fisica","Alcool")) |>
    column_to_rownames("SampleID")
}

socio_r <- rbind(
  load_socio_r("DADOS SOCIODEMOGRÁFICOS - CONTROLE.xlsx", "Controle"),
  load_socio_r("DADOS SOCIODEMOGRÁFICOS - CCL.xlsx",      "CCL")
)

# Variáveis clínicas para correlação
clin_r <- cbind(t(ace_df), socio_r[colnames(t(ace_df)), -1])
clin_r <- apply(clin_r, 2, as.numeric)

# Função: heatmap de correlação Spearman
spearman_heatmap <- function(otu_sub, clin_sub, title, filename) {
  samples_ok <- intersect(colnames(otu_sub), rownames(clin_sub))
  otu_sub    <- otu_sub[, samples_ok]
  clin_sub   <- clin_sub[samples_ok, ]

  cor_mat  <- matrix(NA, nrow(otu_sub), ncol(clin_sub),
                      dimnames = list(rownames(otu_sub), colnames(clin_sub)))
  pval_mat <- cor_mat

  for (i in seq_len(nrow(otu_sub))) {
    for (j in seq_len(ncol(clin_sub))) {
      x <- as.numeric(otu_sub[i, ])
      y <- as.numeric(clin_sub[, j])
      ok <- !is.na(x) & !is.na(y)
      if (sum(ok) >= 5) {
        res <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
        cor_mat[i, j]  <- res$estimate
        pval_mat[i, j] <- res$p.value
      }
    }
  }

  pdf(file.path(OUTDIR, paste0(filename, ".pdf")), width = 13, height = 11)
  pheatmap(cor_mat,
           color  = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101),
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           cluster_cols = FALSE,
           display_numbers = matrix(
             ifelse(pval_mat < 0.001, "***",
             ifelse(pval_mat < 0.01,  "**",
             ifelse(pval_mat < 0.05,  "*",  ""))), nrow = nrow(pval_mat)),
           fontsize_number = 9,
           number_color = "black",
           border_color = NA,
           fontsize_row = 8, fontsize_col = 10,
           main = title)
  dev.off()
  invisible(list(r = cor_mat, p = pval_mat))
}

# Heatmaps
spearman_heatmap(otu_matrix, clin_r,
  sprintf("Spearman OTUs × Variáveis Clínicas — Todos (n=%d)", ncol(otu_matrix)),
  "Spearman_R_All")

spearman_heatmap(otu_matrix[, CTRL_IDS[CTRL_IDS %in% colnames(otu_matrix)]], clin_r,
  sprintf("Spearman OTUs × Variáveis Clínicas — Controle (n=%d)",
          sum(CTRL_IDS %in% colnames(otu_matrix))),
  "Spearman_R_Controle")

spearman_heatmap(otu_matrix[, CCL_IDS[CCL_IDS %in% colnames(otu_matrix)]], clin_r,
  sprintf("Spearman OTUs × Variáveis Clínicas — CCL (n=%d)",
          sum(CCL_IDS %in% colnames(otu_matrix))),
  "Spearman_R_CCL")

message("  Spearman_R_All/Controle/CCL.pdf")

# ============================================================
# PASSO 7: ABUNDÂNCIA DIFERENCIAL — DESeq2 (RECOMENDADO)
# ============================================================
# Requer: BiocManager::install("DESeq2")
# Descomente para usar:

# library(DESeq2)
# ps_counts <- ps   # usar objeto com counts RAW (não rarefado, não relativo)
# dds <- phyloseq_to_deseq2(ps_counts, ~ Grupo)
# dds <- DESeq(dds, test = "Wald", fitType = "parametric")
# res <- results(dds, contrast = c("Grupo","CCL","Controle"))
# res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
#
# # Volcano plot
# res_df <- as.data.frame(res) |>
#   mutate(Significante = padj < 0.05 & abs(log2FoldChange) > 1)
#
# ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significante)) +
#   geom_point(alpha = 0.7, size = 2) +
#   scale_color_manual(values = c("grey70","#C0392B")) +
#   geom_vline(xintercept = c(-1,1), lty = 2, color = "grey50") +
#   geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey50") +
#   labs(title = "Volcano Plot DESeq2 — CCL vs Controle",
#        x = "log2 Fold Change", y = "-log10(p-valor)") +
#   theme_classic()

# ============================================================
# PASSO 8: LEfSe (microbiomeMarker) — OPCIONAL
# ============================================================
# Requer: BiocManager::install("microbiomeMarker")
# Descomente para usar:

# library(microbiomeMarker)
# mm_lefse <- run_lefse(ps, group = "Grupo",
#                       norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)
# plot_ef_bar(mm_lefse) + ggtitle("LEfSe — CCL vs Controle")

# ============================================================
# PASSO 9: CURVAS DE RAREFAÇÃO (vegan)
# ============================================================
message("\n== PASSO 9: Curvas de rarefação ==")

otu_raw <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_raw <- t(otu_raw)

pdf(file.path(OUTDIR, "Rarefaction_Curves_R.pdf"), width = 10, height = 6)
rarecurve(t(otu_raw),
          step    = 1000,
          sample  = min_depth,
          col     = ifelse(colnames(otu_raw) %in% CCL_IDS, "#E07070", "#70B8D4"),
          lwd     = 1.5,
          label   = FALSE,
          main    = "Curvas de Rarefação — CCL (vermelho) vs Controle (azul)",
          xlab    = "Número de reads",
          ylab    = "OTUs observados")
abline(v = min_depth, lty = 2, col = "black", lwd = 1.5)
legend("bottomright",
       legend = c(sprintf("CCL (n=%d)", sum(CCL_IDS %in% colnames(otu_raw))),
                  sprintf("Controle (n=%d)", sum(CTRL_IDS %in% colnames(otu_raw))),
                  sprintf("Rarefação (%s reads)", format(min_depth, big.mark=","))),
       col = c("#E07070","#70B8D4","black"),
       lty = c(1,1,2), lwd = 1.5, cex = 0.9)
dev.off()

message("  Rarefaction_Curves_R.pdf")

# ============================================================
# RESUMO FINAL
# ============================================================
message("\n", paste(rep("=",60), collapse=""))
message("ANÁLISE R CONCLUÍDA")
message(paste(rep("=",60), collapse=""))
message(sprintf("Outputs em: %s", OUTDIR))
message("Figuras:")
message("  AlphaDiversity_Boxplot.pdf/png")
message("  BetaDiversity_PCoA.pdf/png")
message("  Heatmap_R_TopOTUs.pdf/png")
message("  Spearman_R_All/Controle/CCL.pdf")
message("  Rarefaction_Curves_R.pdf")
message("Tabelas:")
message("  Tabela_AlphaDiversity_R.csv")
