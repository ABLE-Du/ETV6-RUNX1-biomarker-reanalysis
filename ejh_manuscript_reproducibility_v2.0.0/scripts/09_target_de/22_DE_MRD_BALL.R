#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# REVISION_01 item 4 -- Differential expression at diagnosis, Poor vs Good D29
# MRD, computed on P_BALL (confirmed B-lineage ALL) instead of the obsolete
# n=455 mixed-lineage cohort.
#
# Cohort : P_BALL, n=191  (see revision_01/reports/LINEAGE_AUDIT.md)
# Counts : integer raw counts, STAR-Counts 'unstranded'
# Design : primary            ~ MRD_group
#          sens_sample_type   ~ sample_type + MRD_group      (BM 134 / PB 57 --
#                              sample type genuinely differs, so this is a
#                              required, not optional, sensitivity model)
#          sens_protocol      ~ protocol + MRD_group          (full cohort)
#          sens_protocol_main ~ protocol + MRD_group          (AALL0232+AALL0331
#                              only; the two 2-patient protocols cause
#                              coefficient separation, same trap as Stage-2)
#          sens_clinical      ~ age_years + sex + log_wbc + MRD_group
#          sens_full          ~ sample_type + protocol + MRD_group (main only)
# Test   : Wald, BH FDR < 0.05, no |log2FC| gate (analysis_lock.yaml)
#
# SHE is NOT modelled here: it fails the locked expression filter (see item 8).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

set.seed(20260915)

ROOT <- Sys.getenv("R_SCRATCH", unset = file.path(tempdir(), "etv6work_r1"))
if (!nzchar(ROOT)) ROOT <- file.path(tempdir(), "etv6work_r1")
ROOT  <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
INTER <- file.path(ROOT, "in")
OUT   <- file.path(ROOT, "out")
RES   <- file.path(OUT, "results", "deseq2")
TAB   <- file.path(OUT, "tables")
FIG   <- file.path(OUT, "figures")
for (d in c(RES, TAB, FIG, file.path(ROOT, "logs"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_lines <- character()
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

w <- max(1, min(6, parallel::detectCores() - 2))
ok <- tryCatch({
  BiocParallel::register(BiocParallel::SnowParam(workers = w))
  TRUE
}, error = function(e) {
  BiocParallel::register(BiocParallel::SerialParam())
  FALSE
})
say(sprintf("[par] workers=%d snow=%s", w, ok))

suppressPackageStartupMessages(library(data.table))

cnt <- data.table::fread(file.path(INTER, "counts_P_BALL.tsv"),
                         sep = "\t", header = TRUE, data.table = FALSE,
                         showProgress = FALSE, colClasses = list(character = 1))
cnt <- as.data.frame(cnt, stringsAsFactors = FALSE)
rownames(cnt) <- cnt[[1]]
cnt <- as.matrix(cnt[, -1, drop = FALSE])
storage.mode(cnt) <- "integer"

meta <- as.data.frame(data.table::fread(file.path(INTER, "meta_P_BALL.tsv"),
                                        sep = "\t", header = TRUE, data.table = FALSE,
                                        showProgress = FALSE, colClasses = "character"),
                      stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
meta <- meta[colnames(cnt), , drop = FALSE]

ann <- as.data.frame(data.table::fread(file.path(INTER, "gene_annotation_aligned.tsv"),
                                       sep = "\t", header = TRUE, data.table = FALSE,
                                       showProgress = FALSE, colClasses = "character"),
                     stringsAsFactors = FALSE)

say(sprintf("[load] P_BALL counts %d genes x %d samples", nrow(cnt), ncol(cnt)))
say(sprintf("[load] MRD_group: %s", paste(names(table(meta$MRD_group)),
                                          table(meta$MRD_group), sep = "=", collapse = ", ")))
say(sprintf("[load] sample_type: %s", paste(names(table(meta$expr_sample_type)),
                                            table(meta$expr_sample_type), sep = "=", collapse = ", ")))
say(sprintf("[load] protocol: %s", paste(names(table(meta$protocol)),
                                         table(meta$protocol), sep = "=", collapse = ", ")))

# ------------------------------------------------- locked pre-filter --------
MIN_COUNT <- 10
MIN_FRAC  <- 0.80
n_need <- ceiling(MIN_FRAC * ncol(cnt))
keep <- rowSums(cnt >= MIN_COUNT) >= n_need
say(sprintf("[filter] counts>=%d in >=%d samples -> %d of %d genes retained (B-ALL universe)",
            MIN_COUNT, n_need, sum(keep), nrow(cnt)))

# ------------------------------------------------------- meta preparation ---
prep <- function(m) {
  m$MRD_group <- factor(m$MRD_group, levels = c("Good", "Poor"))   # Good = reference
  m$sample_type_s <- factor(m$expr_sample_type)
  m$protocol <- factor(m$protocol)
  m$sex <- factor(m$sex)
  m$age_years <- as.numeric(m$age_years)
  m$log_wbc <- as.numeric(m$log_wbc)
  m$cns <- factor(m$cns_status_at_diagnosis)
  m
}
meta <- prep(meta)

MAIN_PROTO <- c("AALL0232", "AALL0331")
main_ids <- meta$sample_id[meta$protocol %in% MAIN_PROTO]
say(sprintf("[subset] main protocols (%s): %d samples, %d events(Poor)",
            paste(MAIN_PROTO, collapse = "+"), length(main_ids),
            sum(meta[main_ids, "MRD_group"] == "Poor")))

fit_one <- function(design_f, label, m, c_mat) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = c_mat, colData = m, design = design_f)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds)
  res <- DESeq2::results(dds, contrast = c("MRD_group", "Poor", "Good"),
                         alpha = 0.05, independentFiltering = TRUE)
  df <- as.data.frame(res)
  df$ensembl_gene_id <- rownames(df)
  rownames(df) <- NULL
  df <- df[, c("ensembl_gene_id", "baseMean", "log2FoldChange", "lfcSE",
               "stat", "pvalue", "padj")]
  list(dds = dds, res = df, label = label, n = ncol(c_mat))
}

designs <- list(
  list(f = ~ MRD_group,                             label = "primary",            main = FALSE),
  list(f = ~ sample_type_s + MRD_group,             label = "sens_sample_type",   main = FALSE),
  list(f = ~ protocol + MRD_group,                  label = "sens_protocol",      main = FALSE),
  list(f = ~ protocol + MRD_group,                  label = "sens_protocol_main", main = TRUE),
  list(f = ~ age_years + sex + log_wbc + MRD_group, label = "sens_clinical",      main = FALSE),
  list(f = ~ sample_type_s + protocol + MRD_group,  label = "sens_full",          main = TRUE)
)

out <- list()
for (d in designs) {
  t0 <- Sys.time()
  if (d$main) {
    mm <- meta[main_ids, , drop = FALSE]; mm$protocol <- droplevels(mm$protocol)
    cm <- cnt[keep, main_ids, drop = FALSE]
  } else {
    mm <- meta; cm <- cnt[keep, , drop = FALSE]
  }
  r <- fit_one(d$f, d$label, mm, cm)
  r$res <- r$res %>% dplyr::left_join(
      ann %>% dplyr::select(ensembl_gene_id, gene_symbol, gene_type), by = "ensembl_gene_id") %>%
    dplyr::select(ensembl_gene_id, gene_symbol, gene_type, dplyr::everything())
  out[[d$label]] <- r
  n_sig <- sum(r$res$padj < 0.05, na.rm = TRUE)
  say(sprintf("[fit] %-20s n=%3d  %5d tested, %4d FDR<0.05, %.1fs",
              d$label, r$n, sum(!is.na(r$res$padj)), n_sig,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

primary <- out$primary
dds <- primary$dds

# ------------------------------------------------------------ write DEGs ----
res_all <- primary$res[order(primary$res$padj, primary$res$pvalue, na.last = TRUE), ]
readr::write_tsv(res_all, file.path(RES, "DEG_MRD_ALL_GENES.tsv.gz"))
say(sprintf("[write] DEG_MRD_ALL_GENES.tsv.gz (%d rows)", nrow(res_all)))

sig <- res_all %>% dplyr::filter(padj < 0.05)
readr::write_tsv(sig, file.path(RES, "DEG_MRD_FDR05.tsv.gz"))
say(sprintf("[write] FDR<0.05: %d genes (%d up in Poor, %d down)",
            nrow(sig), sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))

for (nm in names(out)) {
  x <- out[[nm]]$res[order(out[[nm]]$res$padj, out[[nm]]$res$pvalue, na.last = TRUE), ]
  readr::write_tsv(x, file.path(RES, sprintf("DEG_MRD_%s_ALL_GENES.tsv.gz", nm)))
}
say("[write] per-design DEG tables")

conc <- do.call(rbind, lapply(names(out), function(nm) {
  d <- out[[nm]]$res
  data.frame(model = nm, n_samples = out[[nm]]$n,
             n_tested = sum(!is.na(d$padj)),
             n_FDR05 = sum(d$padj < 0.05, na.rm = TRUE),
             n_up = sum(d$padj < 0.05 & d$log2FoldChange > 0, na.rm = TRUE),
             n_down = sum(d$padj < 0.05 & d$log2FoldChange < 0, na.rm = TRUE),
             median_abs_lfc = round(median(abs(d$log2FoldChange), na.rm = TRUE), 4),
             stringsAsFactors = FALSE)
}))
readr::write_tsv(conc, file.path(TAB, "DE_model_comparison_BALL.tsv"))
say("[write] DE_model_comparison_BALL.tsv")
print(conc)

# ------------------------------------------------- old vs new concordance ---
# The Stage-2 n=455 result is obsolete for B-ALL inference but is retained so we
# can quantify how far it drifted.  Nothing from it is reused as a statistic.
old_path <- file.path(INTER, "DEG_MRD_ALL_GENES_stage2_n455.tsv")
if (file.exists(old_path)) {
  old <- as.data.frame(data.table::fread(old_path, sep = "\t", header = TRUE,
                                         data.table = FALSE, showProgress = FALSE,
                                         colClasses = list(character = 1)))
  new <- res_all
  mrg <- dplyr::inner_join(
    old %>% dplyr::select(ensembl_gene_id, old_lfc = log2FoldChange,
                          old_p = pvalue, old_padj = padj),
    new %>% dplyr::select(ensembl_gene_id, new_lfc = log2FoldChange,
                          new_p = pvalue, new_padj = padj),
    by = "ensembl_gene_id")
  cmp <- data.frame(
    metric = c("n_genes_compared",
               "spearman_lfc", "pearson_lfc",
               "sign_agreement_lfc",
               "spearman_wald_stat",
               "old_FDR05", "new_FDR05", "overlap_FDR05",
               "jaccard_FDR05",
               "spearman_lfc_FDR05_in_new"),
    value = c(
      nrow(mrg),
      round(cor(mrg$old_lfc, mrg$new_lfc, method = "spearman", use = "complete.obs"), 4),
      round(cor(mrg$old_lfc, mrg$new_lfc, method = "pearson", use = "complete.obs"), 4),
      round(mean(sign(mrg$old_lfc) == sign(mrg$new_lfc), na.rm = TRUE), 4),
      NA_real_,
      sum(old$padj < 0.05, na.rm = TRUE),
      sum(new$padj < 0.05, na.rm = TRUE),
      sum(mrg$old_padj < 0.05 & mrg$new_padj < 0.05, na.rm = TRUE),
      round(sum(mrg$old_padj < 0.05 & mrg$new_padj < 0.05, na.rm = TRUE) /
              max(1, sum(mrg$old_padj < 0.05 | mrg$new_padj < 0.05, na.rm = TRUE)), 4),
      round(cor(mrg$old_lfc[mrg$new_padj < 0.05], mrg$new_lfc[mrg$new_padj < 0.05],
                method = "spearman", use = "complete.obs"), 4)))
  readr::write_tsv(cmp, file.path(TAB, "old_vs_BALL_DE_concordance.tsv"))
  say("[write] old_vs_BALL_DE_concordance.tsv")
  print(cmp)
} else {
  say("[warn] Stage-2 n=455 DEG table not pushed; skipping old-vs-new concordance")
}

# ------------------------------------------------------------- vst ---------
vsd <- DESeq2::vst(dds, blind = FALSE)
vst_mat <- as.data.frame(SummarizedExperiment::assay(vsd))
vst_mat$ensembl_gene_id <- rownames(vst_mat)
readr::write_tsv(vst_mat, file.path(INTER, "vst_P_BALL.tsv.gz"))
say(sprintf("[write] vst_P_BALL.tsv.gz %dx%d", nrow(vst_mat) - 1, ncol(vst_mat) - 1))

# ---------------------------------------------------------- QC figures -----
pca <- DESeq2::plotPCA(vsd, intgroup = "MRD_group", returnData = TRUE)
pca$MRD_group <- meta[colnames(vsd), "MRD_group"]
pca$protocol <- meta[colnames(vsd), "protocol"]
p1 <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, colour = MRD_group, shape = protocol)) +
  ggplot2::geom_point(size = 1.8, alpha = 0.8) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "PCA (vst) -- P_BALL (confirmed B-ALL)",
                subtitle = sprintf("n = %d  (Good %d / Poor %d)", ncol(vsd),
                                   sum(meta$MRD_group == "Good"), sum(meta$MRD_group == "Poor")))
ggplot2::ggsave(file.path(FIG, "FigS_PCA_MRD_BALL.pdf"), p1, width = 6.6, height = 5)

rp <- res_all %>% dplyr::filter(!is.na(padj))
p2 <- ggplot2::ggplot(rp, ggplot2::aes(log2FoldChange, -log10(padj))) +
  ggplot2::geom_point(size = 0.5, alpha = 0.35, colour = "grey40") +
  ggplot2::geom_point(data = dplyr::filter(rp, padj < 0.05), size = 0.7,
                      colour = "#B2182B", alpha = 0.8) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "P_BALL diagnosis: Poor vs Good D29 MRD",
                subtitle = sprintf("n=%d, %d genes tested, %d FDR<0.05",
                                   ncol(cnt), nrow(rp), sum(rp$padj < 0.05)))
ggplot2::ggsave(file.path(FIG, "Fig2_volcano_MRD_BALL.pdf"), p2, width = 6.2, height = 5)
say("[write] QC figures")

# ---------------------------------------------------------- provenance -----
session <- utils::sessionInfo()
pv <- data.frame(
  item = c("script", "R_version", "DESeq2_version", "cohort", "n_samples",
           "n_genes_tested", "pre_filter", "seed", "timestamp",
           "lineage_gate", "obsolete_stage2_reused"),
  value = c("analysis/scripts/22_DE_MRD_BALL.R", session$R.version$version.string,
            as.character(packageVersion("DESeq2")), "P_BALL (confirmed B-lineage)", ncol(cnt),
            nrow(res_all), sprintf("counts>=%d in >=%.0f%% samples", MIN_COUNT, 100 * MIN_FRAC),
            "20260915", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
            "lineage == B-ALL (patient-level annotation)", "no"))
readr::write_tsv(pv, file.path(RES, "DE_MRD_BALL_provenance.tsv"))

writeLines(log_lines, file.path(ROOT, "logs", "22_DE_MRD_BALL.log"))
say("[done]")
