############################################################
# RRBS methylation analysis (methylKit) — GitHub-ready template
# Human astrocytes; ATP13A2 LOF (c.1306) ± CGP
#
# Covers:
#  - Import methylation calls (bedGraph or Bismark coverage)
#  - QC: methylation/coverage distributions
#  - Filtering by coverage
#  - Unite CpGs across samples
#  - Correlation, clustering, PCA
#  - Global methylation (descriptive)
#  - Regional methylation via tiling windows
#  - Differential methylation (logistic regression w/ overdispersion)
#  - Extraction of hyper/hypo DMRs (q + diff cutoffs)
#  - bedGraph export for genome browser
#  - UCSC/Table Browser export (BED) for gene annotation workflow
#  - OPTIONAL: integration with expression + enrichment (clusterProfiler)
#
# Author: <>
# Last updated: <>
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(tidyverse)
  library(data.table)
  library(ggplot2)
})

# Optional packages (only used if you run those sections)
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
has_cp       <- requireNamespace("clusterProfiler", quietly = TRUE)
has_orgdb    <- requireNamespace("org.Hs.eg.db", quietly = TRUE)

set.seed(1)

# ==========================================================
# 0) CONFIG (edit these)
# ==========================================================

# Project layout
PROJECT_DIR <- getwd()
META_PATH   <- file.path(PROJECT_DIR, "metadata", "samplesheet.csv")
OUT_DIR     <- file.path(PROJECT_DIR, "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "qc"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "similarity"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "global"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tiles"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "diffmeth"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "browser_tracks"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "exports_ucsc"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "integration"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "enrichment"), showWarnings = FALSE, recursive = TRUE)

# Genome build
GENOME_ASSEMBLY <- "hg38"

# Input format:
# - If you have Bismark coverage (*.cov / *.bismark.cov): use "bismarkCoverage"
# - If you truly have 4-column bedGraph (chr, start, end, %meth): methylKit cannot
#   compute logistic regression without counts, so you MUST convert bedGraph -> 6-col
#   coverage (or provide counts another way). Conversion helper included below.
PIPELINE_TYPE <- "bismarkCoverage"  # change if needed

# Coverage filtering thresholds
MIN_COV     <- 10      # per-sample minimum coverage per CpG (recommended in methylKit paper)
MAX_COV_PCT <- 99.9    # remove extreme high coverage (PCR bias / clonal reads)

# Unite rules (CpGs shared across samples)
# NOTE: methylKit has min.per.group but it assumes a treatment vector. We’ll handle
# pairwise contrasts by subsetting. unite() keeps sites across all samples by default.
UNITE_MINCOV <- MIN_COV

# Tiling params (regional DMRs)
TILE_SIZE_BP <- 1000
TILE_STEP_BP <- 1000   # contiguous windows if equal to TILE_SIZE_BP
TILE_COVBASES <- MIN_COV

# Differential methylation modeling
OVERDISPERSION_MODE <- "MN"   # "none", "MN", "shrinkMN"
ADJ_METHOD          <- "BH"   # "BH" (FDR) or "SLIM" (methylKit SLIM)

# Default thresholds for DMR extraction (from your methods)
Q_CUTOFF <- 0.01
DIFF_STRICT <- 75
DIFF_RELAX  <- 25

# ==========================================================
# 1) OPTIONAL HELPER: convert 4-col bedGraph -> 6-col coverage
# ==========================================================
# If your bedGraph only has % methylation and NO counts, you cannot do logistic
# regression (needs methylated/unmethylated counts).
#
# If you *also* have total coverage per site in some file, you can reconstruct counts:
#   meth_count = round(coverage * pct/100)
#   unmeth_count = coverage - meth_count
#
# This helper expects:
#   bedgraph: chr, start, end, pct_meth
#   coverage: chr, start, end, coverage
# and writes a 6-col Bismark-like coverage file:
#   chr, start, end, pct_meth, meth_count, unmeth_count
convert_bedgraph_pct_to_cov6 <- function(bedgraph_pct_path,
                                        coverage_path,
                                        out_cov6_path) {
  bg <- fread(bedgraph_pct_path, header = FALSE)
  cov <- fread(coverage_path, header = FALSE)

  stopifnot(ncol(bg) >= 4, ncol(cov) >= 4)
  setnames(bg, 1:4, c("chr", "start", "end", "pct"))
  setnames(cov, 1:4, c("chr", "start", "end", "coverage"))

  m <- merge(bg, cov, by = c("chr", "start", "end"))
  m[, meth := as.integer(round(coverage * pct / 100))]
  m[, unmeth := as.integer(pmax(coverage - meth, 0))]
  m[, pct := as.numeric(pct)]

  out <- m[, .(chr, start, end, pct, meth, unmeth)]
  fwrite(out, out_cov6_path, sep = "\t", col.names = FALSE)
  return(invisible(out_cov6_path))
}

# ==========================================================
# 2) Load samplesheet
# ==========================================================
# samplesheet.csv REQUIRED columns:
#   sample_id, condition, file_path
#
# Example conditions: WT, c1306, c1306_CGP
samples <- readr::read_csv(META_PATH, show_col_types = FALSE) %>%
  mutate(file_path = ifelse(file.exists(file_path),
                            file_path,
                            file.path(PROJECT_DIR, file_path)))

stopifnot(all(c("sample_id", "condition", "file_path") %in% colnames(samples)))
stopifnot(all(file.exists(samples$file_path)))

readr::write_csv(samples, file.path(OUT_DIR, "samplesheet_resolved.csv"))

# Dummy treatment at import (we set treatments per-contrast later)
dummy_treatment <- rep(0, nrow(samples))

# ==========================================================
# 3) Import methylation calls (CpG context)
# ==========================================================
message("Reading methylation calls with pipeline = ", PIPELINE_TYPE)

meth_raw <- methRead(
  location  = samples$file_path,
  sample.id = samples$sample_id,
  assembly  = GENOME_ASSEMBLY,
  treatment = dummy_treatment,
  context   = "CpG",
  pipeline  = PIPELINE_TYPE,
  mincov    = 1
)

# ==========================================================
# 4) QC plots: coverage + methylation distributions
# ==========================================================
qc_dir <- file.path(OUT_DIR, "qc")

for (i in seq_along(meth_raw)) {
  sid <- samples$sample_id[i]

  pdf(file.path(qc_dir, paste0("coverage_stats_", sid, ".pdf")), width = 7, height = 5)
  getCoverageStats(meth_raw[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()

  pdf(file.path(qc_dir, paste0("methylation_stats_", sid, ".pdf")), width = 7, height = 5)
  getMethylationStats(meth_raw[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
}

qc_table <- tibble(
  sample_id = samples$sample_id,
  condition = samples$condition,
  n_CpGs_imported = sapply(meth_raw, nrow),
  mean_coverage   = sapply(meth_raw, function(x) mean(getData(x)$coverage, na.rm = TRUE)),
  median_coverage = sapply(meth_raw, function(x) median(getData(x)$coverage, na.rm = TRUE))
)

readr::write_csv(qc_table, file.path(qc_dir, "qc_table_basic_methylkit.csv"))

# ==========================================================
# 5) Filter by coverage (remove low + extreme high)
# ==========================================================
meth_filt <- filterByCoverage(
  meth_raw,
  lo.count = MIN_COV,
  hi.perc  = MAX_COV_PCT
)

saveRDS(meth_filt, file.path(OUT_DIR, "meth_filtered_MethylRawList.rds"))

# ==========================================================
# 6) Unite CpGs across samples
# ==========================================================
meth_base <- unite(
  meth_filt,
  destrand = FALSE,
  mincov   = UNITE_MINCOV
)

saveRDS(meth_base, file.path(OUT_DIR, "meth_united_CpG_MethylBase.rds"))

# ==========================================================
# 7) Sample similarity: correlation, clustering, PCA
# ==========================================================
sim_dir <- file.path(OUT_DIR, "similarity")

# Correlation matrix (only sites covered in all samples internally)
cor_mat <- getCorrelation(meth_base, plot = FALSE, method = "pearson")
write.csv(cor_mat, file.path(sim_dir, "correlation_matrix_pearson.csv"))

png(file.path(sim_dir, "correlation_heatmap_pearson.png"), width = 900, height = 800, res = 150)
if (has_pheatmap) {
  pheatmap::pheatmap(cor_mat)
} else {
  heatmap(cor_mat)
}
dev.off()

# Hierarchical clustering (distance = 1 - correlation recommended in paper)
png(file.path(sim_dir, "cluster_samples.png"), width = 900, height = 700, res = 150)
clusterSamples(meth_base,
               dist = "correlation",
               method = "ward.D2",
               plot = TRUE)
dev.off()

# PCA
pca_obj <- PCASamples(meth_base, obj.return = TRUE)

pca_df <- as_tibble(pca_obj$x[, 1:3], rownames = "sample_id") %>%
  left_join(samples %>% select(sample_id, condition), by = "sample_id")

readr::write_csv(pca_df, file.path(sim_dir, "pca_scores_PC1_PC3.csv"))

p_pca12 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample_id)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(max.overlaps = 50) +
  theme_classic() +
  labs(title = "PCA of methylation profiles (sites shared across samples)")
ggsave(file.path(sim_dir, "PCA_PC1_PC2.png"), p_pca12, width = 7, height = 5, dpi = 300)

# ==========================================================
# 8) Global methylation (descriptive)
# ==========================================================
pm <- percMethylation(meth_base)              # CpGs x samples matrix
global_pct <- colMeans(pm, na.rm = TRUE)

global_df <- tibble(
  sample_id = colnames(pm),
  global_pct_methylation = as.numeric(global_pct)
) %>%
  left_join(samples %>% select(sample_id, condition), by = "sample_id")

readr::write_csv(global_df, file.path(OUT_DIR, "global", "global_methylation_by_sample.csv"))

p_global <- ggplot(global_df, aes(x = condition, y = global_pct_methylation)) +
  geom_violin(trim = FALSE) +
  geom_point(size = 2, position = position_jitter(width = 0.1)) +
  theme_classic() +
  labs(x = NULL, y = "Global CpG methylation (%)")
ggsave(file.path(OUT_DIR, "global", "global_methylation_violin.png"),
       p_global, width = 6, height = 4, dpi = 300)

# ==========================================================
# 9) Regional methylation: tiling windows
# ==========================================================
tiles_raw <- tileMethylCounts(
  meth_filt,
  win.size  = TILE_SIZE_BP,
  step.size = TILE_STEP_BP,
  cov.bases = TILE_COVBASES
)

tiles_base <- unite(
  tiles_raw,
  destrand = FALSE,
  mincov   = MIN_COV
)

saveRDS(tiles_base, file.path(OUT_DIR, "tiles", "meth_united_TILES_MethylBase.rds"))

# ==========================================================
# 10) Differential methylation (pairwise contrasts)
# ==========================================================
# methylKit chooses logistic regression when there are replicates;
# Fisher's exact test is used if one sample per group.
#
# We run logistic regression via calculateDiffMeth(), with overdispersion correction.
#
# Note: treatment vector defines control vs treatment:
#  - treatment = 0 => control/reference
#  - treatment = 1 => experimental/test

run_pairwise_diffmeth <- function(meth_obj, samples_df, control_cond, test_cond,
                                  out_prefix,
                                  q_cutoff = Q_CUTOFF,
                                  diff_cutoff = DIFF_RELAX,
                                  overdispersion = OVERDISPERSION_MODE,
                                  adjust_method = ADJ_METHOD,
                                  n_cores = 1) {

  outdir <- file.path(OUT_DIR, "diffmeth", out_prefix)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  keep <- samples_df$condition %in% c(control_cond, test_cond)
  s2 <- samples_df[keep, , drop = FALSE]

  # Ensure methylKit object is reorganized to match subset sample order
  meth_sub <- reorganize(meth_obj, sample.ids = s2$sample_id)

  tr <- ifelse(s2$condition == test_cond, 1, 0)
  treatment(meth_sub) <- tr

  diff_obj <- calculateDiffMeth(
    meth_sub,
    overdispersion = overdispersion,
    adjust         = adjust_method,
    mc.cores       = n_cores
  )

  diff_tbl <- as_tibble(getData(diff_obj))

  # Extract significant calls (hypermethylated/hypomethylated relative to control)
  sig_tbl <- getMethylDiff(
    diff_obj,
    difference = diff_cutoff,
    qvalue     = q_cutoff
  ) %>%
    getData() %>%
    as_tibble() %>%
    mutate(direction = ifelse(meth.diff > 0, "hypermethylated", "hypomethylated")) %>%
    arrange(qvalue)

  readr::write_csv(diff_tbl, file.path(outdir, paste0(out_prefix, "_ALL.csv")))
  readr::write_csv(sig_tbl,  file.path(outdir, paste0(out_prefix, "_SIG_q", q_cutoff, "_diff", diff_cutoff, ".csv")))

  # Volcano-like plot: meth.diff vs -log10(q)
  p_vol <- ggplot(diff_tbl, aes(x = meth.diff, y = -log10(qvalue))) +
    geom_point(alpha = 0.25) +
    geom_hline(yintercept = -log10(q_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-diff_cutoff, diff_cutoff), linetype = "dashed") +
    theme_classic() +
    labs(
      title = out_prefix,
      x = "Methylation difference (%), test - control",
      y = "-log10(q-value)"
    )
  ggsave(file.path(outdir, paste0(out_prefix, "_volcano.png")), p_vol, width = 6, height = 5, dpi = 300)

  # Chromosome-level summary: % hyper vs hypo among significant calls
  chr_summary <- sig_tbl %>%
    count(chr, direction) %>%
    group_by(chr) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  readr::write_csv(chr_summary, file.path(outdir, paste0(out_prefix, "_chr_hyper_hypo_pct.csv")))

  p_chr <- ggplot(chr_summary, aes(x = chr, y = pct, fill = direction)) +
    geom_col(position = "stack") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0(out_prefix, ": hyper/hypo composition by chromosome"),
         x = NULL, y = "Percent of significant calls")

  ggsave(file.path(outdir, paste0(out_prefix, "_chr_hyper_hypo_pct.png")),
         p_chr, width = 10, height = 4, dpi = 300)

  # Export BED for UCSC Table Browser annotation (recommended manual step)
  bed <- sig_tbl %>%
    transmute(
      chr   = chr,
      start = start,
      end   = end,
      name  = paste0(out_prefix, "_", row_number()),
      score = 0,
      strand = if ("strand" %in% colnames(sig_tbl)) strand else "*"
    )

  bed_path <- file.path(OUT_DIR, "exports_ucsc", paste0(out_prefix, "_SIG.bed"))
  data.table::fwrite(bed, bed_path, sep = "\t", col.names = FALSE)

  # Export bedGraph track for browser visualization
  # methylKit expects diff_obj and will output differential methylation track
  # (this is a differential track, not raw methylation per sample)
  bedgraph_path <- file.path(OUT_DIR, "browser_tracks", paste0(out_prefix, "_diffmeth.bedGraph"))
  bedgraph(diff_obj, file.name = bedgraph_path, col.name = out_prefix)

  return(list(diff_obj = diff_obj, diff_tbl = diff_tbl, sig_tbl = sig_tbl,
              bed_path = bed_path, bedgraph_path = bedgraph_path))
}

# ---- CpG-level contrasts ----
res_cpg_c1306_vs_WT_strict <- run_pairwise_diffmeth(
  meth_obj      = meth_base,
  samples_df    = samples,
  control_cond  = "WT",
  test_cond     = "c1306",
  out_prefix    = "CpG_c1306_vs_WT_STRICT",
  q_cutoff      = Q_CUTOFF,
  diff_cutoff   = DIFF_STRICT
)

res_cpg_CGP_vs_c1306_relax <- run_pairwise_diffmeth(
  meth_obj      = meth_base,
  samples_df    = samples,
  control_cond  = "c1306",
  test_cond     = "c1306_CGP",
  out_prefix    = "CpG_c1306_CGP_vs_c1306_RELAX",
  q_cutoff      = Q_CUTOFF,
  diff_cutoff   = DIFF_RELAX
)

# ---- TILE-level contrasts (regional DMRs) ----
res_tile_c1306_vs_WT_strict <- run_pairwise_diffmeth(
  meth_obj      = tiles_base,
  samples_df    = samples,
  control_cond  = "WT",
  test_cond     = "c1306",
  out_prefix    = "TILE_c1306_vs_WT_STRICT",
  q_cutoff      = Q_CUTOFF,
  diff_cutoff   = DIFF_STRICT
)

res_tile_CGP_vs_c1306_relax <- run_pairwise_diffmeth(
  meth_obj      = tiles_base,
  samples_df    = samples,
  control_cond  = "c1306",
  test_cond     = "c1306_CGP",
  out_prefix    = "TILE_c1306_CGP_vs_c1306_RELAX",
  q_cutoff      = Q_CUTOFF,
  diff_cutoff   = DIFF_RELAX
)

# ==========================================================
# 11) OPTIONAL: integrate methylation with expression (requires annotated genes)
# ==========================================================
# Workflow:
# 1) Upload BED from results/exports_ucsc/ to UCSC Table Browser (hg38, knownCanonical)
# 2) Download annotation table with gene_symbol column
# 3) Merge annotation back onto methylKit significant results by chr/start/end
# 4) Provide DEG file with gene_symbol + log2FC (or your column names)

integrate_meth_expr <- function(sig_tbl, deg_path,
                                gene_col = "gene_symbol",
                                lfc_col  = "log2FC",
                                out_prefix = "integration") {

  deg <- readr::read_csv(deg_path, show_col_types = FALSE)
  stopifnot(gene_col %in% colnames(sig_tbl))
  stopifnot(gene_col %in% colnames(deg), lfc_col %in% colnames(deg))

  m <- sig_tbl %>%
    inner_join(deg %>% select(all_of(gene_col), all_of(lfc_col)), by = gene_col) %>%
    rename(log2FC = all_of(lfc_col)) %>%
    mutate(direction = ifelse(meth.diff > 0, "hypermethylated", "hypomethylated"))

  outdir <- file.path(OUT_DIR, "integration", out_prefix)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  readr::write_csv(m, file.path(outdir, paste0(out_prefix, "_meth_expr_merged.csv")))

  # Scatterplots hyper/hypo
  p_hyper <- ggplot(filter(m, direction == "hypermethylated"), aes(x = meth.diff, y = log2FC)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(title = paste0(out_prefix, ": hyper"), x = "Methylation difference (%)", y = "Expression log2FC")

  p_hypo <- ggplot(filter(m, direction == "hypomethylated"), aes(x = meth.diff, y = log2FC)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(title = paste0(out_prefix, ": hypo"), x = "Methylation difference (%)", y = "Expression log2FC")

  ggsave(file.path(outdir, paste0(out_prefix, "_scatter_hyper.png")), p_hyper, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outdir, paste0(out_prefix, "_scatter_hypo.png")),  p_hypo,  width = 6, height = 5, dpi = 300)

  # Heatmap: top 10 hyper + top 10 hypo
  if (has_pheatmap) {
    top_hyper <- m %>% filter(direction == "hypermethylated") %>% arrange(desc(meth.diff)) %>% slice_head(n = 10)
    top_hypo  <- m %>% filter(direction == "hypomethylated") %>% arrange(meth.diff) %>% slice_head(n = 10)
    top <- bind_rows(top_hyper, top_hypo) %>% distinct(.data[[gene_col]], .keep_all = TRUE)

    mat <- top %>%
      select(!!gene_col, meth.diff, log2FC) %>%
      rename(gene = !!gene_col) %>%
      column_to_rownames("gene") %>%
      as.matrix()

    png(file.path(outdir, paste0(out_prefix, "_heatmap_top10_hyper_hypo.png")), width = 900, height = 700, res = 150)
    pheatmap::pheatmap(mat, scale = "none")
    dev.off()
  }

  return(invisible(m))
}

# Example (uncomment after you have annotated sig tables with gene symbols):
# integrate_meth_expr(
#   sig_tbl    = res_tile_c1306_vs_WT_strict$sig_tbl,   # must contain gene_symbol
#   deg_path   = "data/expression_DEG.csv",
#   gene_col   = "gene_symbol",
#   lfc_col    = "log2FC",
#   out_prefix = "TILE_c1306_vs_WT"
# )

# ==========================================================
# 12) OPTIONAL: enrichment with clusterProfiler
# ==========================================================
run_enrichment_go <- function(gene_symbols, out_prefix, ont = "BP") {
  if (!has_cp || !has_orgdb) {
    message("clusterProfiler/org.Hs.eg.db not installed; skipping enrichment.")
    return(NULL)
  }
  library(clusterProfiler)
  library(org.Hs.eg.db)

  outdir <- file.path(OUT_DIR, "enrichment", out_prefix)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  eg <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(eg) == 0) return(NULL)

  ego <- enrichGO(
    gene          = eg$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  ego_df <- as.data.frame(ego)
  readr::write_csv(ego_df, file.path(outdir, paste0(out_prefix, "_enrichGO_", ont, ".csv")))

  if (nrow(ego_df) > 0) {
    png(file.path(outdir, paste0(out_prefix, "_dotplot.png")), width = 1100, height = 800, res = 150)
    print(dotplot(ego, showCategory = 20))
    dev.off()
  }
  return(ego)
}

# Example (uncomment after your sig tables have gene_symbol column):
# hyper_genes <- unique(res_tile_c1306_vs_WT_strict$sig_tbl$gene_symbol[res_tile_c1306_vs_WT_strict$sig_tbl$meth.diff > 0])
# run_enrichment_go(hyper_genes, "TILE_c1306_vs_WT_hyper_BP", ont = "BP")

# ==========================================================
# 13) Save session info
# ==========================================================
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("All done. Outputs in: ", OUT_DIR)
