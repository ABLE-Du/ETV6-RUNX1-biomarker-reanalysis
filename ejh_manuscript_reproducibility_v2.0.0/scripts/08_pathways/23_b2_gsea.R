#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 23_b2_gsea.R -- B2 step 1: pre-ranked GSEA on the OKSA genome-scale ranking.
#
# Ranking: the EOI slow-vs-fast coefficient for EVERY gene in Oksa Supplementary
# Table 22 (19,588 rows).  No DEG filtering -- the brief forbids DEG-only GSEA and
# a filtered ranking is exactly what produces the "only significant pathways are
# reported" failure mode.
#
# Collections: Hallmark, Reactome, GO:BP, KEGG (MSigDB via msigdbr).
#
# Duplicate gene symbols in the Oksa table (3,660 rows) are collapsed by the MEAN
# coefficient within a symbol, and the number affected is written to the QC file so
# the decision is visible rather than implicit.
#
# ALL pathways are saved with NES, p, padj and leadingEdge.  Significance is an
# additional column, never a filter.
#
# Outputs (under out/, pulled back into the project)
#   results/B2_oksa_fgsea_all_pathways.tsv.gz
#   intermediate/B2_oksa_fgsea_qc.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(fgsea); library(msigdbr)
})

ROOT <- Sys.getenv("B2_ROOT", unset = ".")
IN <- file.path(ROOT, "intermediate", "oksa")
OUTR <- file.path(ROOT, "results"); OUTI <- file.path(ROOT, "intermediate")
dir.create(OUTR, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(paste0(...), "\n")

oksa <- fread(file.path(IN, "SuppTable22_DE_analyses.tsv"),
              sep = "\t", quote = "", na.strings = character(0))
coef_col <- "EOI slow VS fast | Coefficient slow VS fast"
p_col <- "EOI slow VS fast | P-value"
stopifnot(coef_col %in% names(oksa))
oksa[, symbol := toupper(trimws(get("Gene symbol")))]
oksa[, coef := suppressWarnings(as.numeric(get(coef_col)))]
oksa[, pval := suppressWarnings(as.numeric(get(p_col)))]

n_rows <- nrow(oksa)
n_nocoef <- sum(is.na(oksa$coef))
n_nosym <- sum(oksa$symbol == "" | is.na(oksa$symbol))

d <- oksa[!is.na(oksa$coef) & oksa$symbol != ""]
dup_syms <- unique(d$symbol[duplicated(d$symbol)])
d <- d[, .(coef = mean(coef), n_rows = .N), by = symbol]
say(sprintf("[rank] input rows=%d  usable rows=%d  symbols with a numeric coef=%d",
            n_rows, nrow(d), length(unique(d$symbol))))
say(sprintf("[rank] rows lacking a coefficient=%d  rows lacking a symbol=%d  collapsed duplicate symbols=%d",
            n_nocoef, n_nosym, length(dup_syms)))

ranks <- setNames(d$coef, d$symbol)
ranks <- sort(ranks, decreasing = TRUE)

qc <- data.table(
  item = c("input_rows", "usable_rows", "unique_symbols", "rows_without_coefficient",
           "rows_without_symbol", "collapsed_duplicate_symbols", "ranking_min",
           "ranking_max", "ranking_median"),
  value = c(n_rows, nrow(d), length(unique(d$symbol)), n_nocoef, n_nosym,
            length(dup_syms), min(ranks), max(ranks), median(ranks)))
fwrite(qc, file.path(OUTI, "B2_oksa_fgsea_qc.tsv"), sep = "\t")

# ---------------------------------------------------------------- gene sets
sets <- list()
hall <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) NULL)
if (!is.null(hall)) sets$Hallmark <- hall
sets$Reactome <- msigdbr(species = "Homo sapiens", collection = "C2",
                         subcollection = "CP:REACTOME")
sets$GO_BP <- msigdbr(species = "Homo sapiens", collection = "C5",
                      subcollection = "GO:BP")
sets$KEGG <- msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:KEGG_LEGACY")

all_res <- list()
for (nm in names(sets)) {
  df <- as.data.frame(sets[[nm]])
  gs <- split(toupper(df$gene_symbol), df$gs_name)
  gs <- lapply(gs, function(x) unique(x[!is.na(x) & x != ""]))
  gs <- gs[vapply(gs, length, integer(1)) >= 10 & vapply(gs, length, integer(1)) <= 500]
  say(sprintf("[sets] %-9s pathways kept=%d", nm, length(gs)))
  set.seed(20260918)
  r <- fgsea::fgsea(pathways = gs, stats = ranks, minSize = 10, maxSize = 500,
                    eps = 0, nPermSimple = 10000)
  r$collection <- nm
  r$leadingEdge <- vapply(r$leadingEdge, paste, character(1), collapse = ";")
  all_res[[nm]] <- as.data.table(r)
}

res <- rbindlist(all_res, fill = TRUE)
res[, direction := ifelse(NES > 0, "higher_in_slow", "higher_in_fast")]
res <- res[order(padj, pval)]
fwrite(res, file.path(OUTR, "B2_oksa_fgsea_all_pathways.tsv.gz"), sep = "\t")

say(sprintf("[fgsea] total pathways=%d  padj<0.05=%d  padj<0.25=%d",
            nrow(res), sum(res$padj < 0.05, na.rm = TRUE),
            sum(res$padj < 0.25, na.rm = TRUE)))
say("[done] B2 step 1")
