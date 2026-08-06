#!/usr/bin/env Rscript
# 18_tf_activity.R: COMPUTE. Calibration of the committed TF_Targets sweep's TF ranks.
# =============================================================================
# The top-ranking regulons of the committed TF_Targets sweep are re-scored here against
# four properties of the network they were read from. Sections 6 to 9 take one each:
# the rank's stability across network variant and estimator, which targets carry the
# score and how many other regulons claim them, activity against regulon size with
# matched random-regulon nulls, and whether the recorded per-edge signs change the answer.
# These are the forensics the mouse anchor ran on its own HIF result
# (mouse_anchor/03_results/04_tf), applied to the same question here.
#
# Annotation tier. No row reaches 03_results/master/ or any effect-size accumulator. A
# regulon scoring high says its target genes move with the synovial-fluid side of this
# contrast.
#
# Reads, read-only:
#   03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv     signed moderated t, HGNC symbol
#   03_results/03_pseudobulk/tables/de_SFvsPB_{treg,tcon,cd8}.csv  avg_expr + padj per gene
#   03_results/03_pseudobulk/tables/gene_symbols.csv               the Ensembl-to-HGNC seam map
#   03_results/14_unbiased_enrichment/tables/gsea_all.csv          the committed TF_Targets rows
#   ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv  (SHA-256 pinned)
#
# Writes 03_results/18_tf_activity/tables/, each captioned in that stage's README:
#   source_hash_manifest.csv       ranked_list_keycheck.csv      symbol_vocabulary_probes.csv
#   symbol_vocabulary_check.csv    alias_recovery.csv            network_variants.csv
#   tf_activity_all.csv            fgsea_family_size_cap.csv     hif1a_rank_cascade.csv
#   target_decomposition.csv       target_decomposition_summary.csv
#   regulon_size_calibration.csv   regulon_size_spearman.csv     size_matched_null.csv
#   signed_vs_unsigned.csv         canonical_hif1a_targets.csv
#
# Writes 03_results/objects/ as checkpoints:
#   18_tf_networks.rds   the network variants per population
#   18_tf_activity.rds   raw decoupleR output for every configuration
#
# Compute only. Figures live in 18_tf_activity_viz.R.
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/18_tf_activity.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(decoupleR)
  library(fgsea)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")          # FIG_CFG, round_numeric_cols
source("02_analysis/helpers/source_hash_manifest.R")  # source_sha256, verify_source_hash
source("02_analysis/helpers/symbol_alias.R")          # build_alias_map, ALIAS_RESOLUTIONS

STAGE  <- "18_tf_activity"
SCRIPT <- "02_analysis/scripts/18_tf_activity.R"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# 0. Config: every parameter read from analysis_config.yaml
# ============================================================================

CFG <- FIG_CFG
TA  <- CFG$tf_activity
if (is.null(TA))
  stop("[18] analysis_config.yaml has no `tf_activity:` block — add it before running.")

THR      <- CFG$thresholds
MINSZ    <- as.integer(THR$gsea_min_size %||% 5L)
MAXSZ    <- as.integer(THR$gsea_max_size %||% 500L)
SEED     <- as.integer(THR$gsea_seed     %||% 123L)

NET_PATH <- CFG$unbiased_enrichment$tf_network$path
if (is.null(NET_PATH))
  stop("[18] `unbiased_enrichment.tf_network.path` is unset; the network path is config, never inline.")
MOR_COL  <- TA$mor_col %||% "weight"
SD_COL   <- TA$sign_decision_col %||% "sign_decision"
SD_DEF   <- TA$default_sign_decision %||% "default activation"

VARIANTS <- unlist(TA$network_variants)
METHODS  <- unlist(TA$methods)
CONS_ST  <- unlist(TA$consensus_statistics)
FOCUS    <- unlist(TA$focus_tfs)
DECOMP   <- unlist(TA$decompose_tfs)
PRIMARY  <- TA$primary_population %||% "treg"
SEL_MAX  <- as.integer(TA$selective_max_regulons %||% 1L)
BANDS    <- as.integer(unlist(TA$promiscuity_bands))
NDRAW    <- as.integer(TA$null_draws %||% 2000L)
NDEC     <- as.integer(TA$null_expression_deciles %||% 10L)
CANON    <- unlist(TA$canonical_hif1a_targets)
FG_TABLE <- TA$fgsea_sweep_table
FG_DB    <- TA$fgsea_database %||% "TF_Targets"

POPS <- vapply(TA$populations, function(p) p$label, character(1))
names(POPS) <- vapply(TA$populations, function(p) p$tag, character(1))

RESULTS <- CFG$paths$results %||% "03_results/"
TBL     <- file.path(RESULTS, STAGE, CFG$paths$stage_tables_subdir %||% "tables")
DIR_OBJ <- CFG$paths$objects %||% "03_results/objects/"
for (d in c(TBL, file.path(TBL, "_overview"), DIR_OBJ))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(SEED)

emit <- function(df, fname) {
  p <- file.path(TBL, fname)
  readr::write_csv(round_numeric_cols(df), p)
  message(sprintf("  [SAVE] %-36s %5d rows x %2d cols", fname, nrow(df), ncol(df)))
  invisible(p)
}

message("=================================================================")
message("18_tf_activity — inferred TF activity on the SF-vs-PB niche contrast")
message("=================================================================")

# ============================================================================
# 1. Pin the cross-compartment network
# ============================================================================
# The first run writes the pin and later runs are gated against it, so a network that
# moves under this stage stops the run.

message("\n[1] Pinning the cross-compartment CollecTRI network ...")
MANIFEST  <- file.path(TBL, "source_hash_manifest.csv")
NET_LABEL <- "collectri_regulons_human"
if (!file.exists(NET_PATH))
  stop("[18] CollecTRI network absent: ", NET_PATH)
net_key <- manifest_key(NET_PATH, root = getwd())
if (!file.exists(MANIFEST)) {
  readr::write_csv(tibble(source_label = NET_LABEL, source_path = net_key,
                          sha256 = source_sha256(NET_PATH)), MANIFEST)
  message("  pin written (first run): ", MANIFEST)
}
NET_SHA <- verify_source_hash(NET_PATH, NET_LABEL, MANIFEST, root = getwd())
message(sprintf("  %s  sha256 %s", net_key, NET_SHA))

# ============================================================================
# 2. Ranked lists and the symbol-vocabulary guard
# ============================================================================
# The trap this section guards. The JIA count matrix carries pre-2019 HGNC symbols
# (MB21D1 for CGAS, TMEM173 for STING1, MARCH5 for MARCHF5, MRE11A for MRE11) while
# CollecTRI carries current ones, so a renamed target drops out of a regulon silently.
# Every unmatched target is therefore reported with its cause: absent from the count
# matrix, dropped by the expression filter, or recoverable through an HGNC alias.

message("\n[2] Ranked lists, expression annotation, and the key check ...")

RANK_TPL <- CFG$coresh$query_signatures$ranked %||%
  "03_results/03_pseudobulk/tables/ranked_{population}.tsv"
DE_TPL   <- CFG$coresh$query_signatures$de_table %||%
  "03_results/03_pseudobulk/tables/de_SFvsPB_{population}.csv"
fill_tpl <- function(tpl, tag) gsub("{population}", tag, tpl, fixed = TRUE)

read_ranked <- function(tag) {
  p <- fill_tpl(RANK_TPL, tag)
  stopifnot(file.exists(p))
  d <- readr::read_tsv(p, col_names = c("gene", "stat"), show_col_types = FALSE)
  stopifnot(!any(duplicated(d$gene)))
  setNames(d$stat, d$gene)
}
read_de <- function(tag) {
  p <- fill_tpl(DE_TPL, tag)
  stopifnot(file.exists(p))
  readr::read_csv(p, show_col_types = FALSE)
}

STATS <- lapply(names(POPS), read_ranked); names(STATS) <- names(POPS)
DE    <- lapply(names(POPS), read_de);     names(DE)    <- names(POPS)

# The seam map holds every symbol the count matrix carried before filterByExpr, which is
# what separates a renamed symbol from one dropped for low expression below.
MATRIX_SYMBOLS <- unique(readr::read_csv(
  "03_results/03_pseudobulk/tables/gene_symbols.csv", show_col_types = FALSE)$gene_symbol)

# Ensembl ids in a ranked list intersect every reference at ~zero, and fgsea and decoupleR
# report that as an empty result. Guard before anything runs.
keycheck <- bind_rows(lapply(names(POPS), function(tag) {
  g <- names(STATS[[tag]])
  tibble(population = POPS[[tag]], n_genes = length(g),
         n_ensembl_like = sum(grepl("^ENSG[0-9]{6,}", g)),
         frac_ensembl_like = mean(grepl("^ENSG[0-9]{6,}", g)),
         key = if (mean(grepl("^ENSG[0-9]{6,}", g)) > 0.5) "ensembl_id" else "hgnc_symbol")
}))
print(as.data.frame(keycheck), row.names = FALSE)
if (any(keycheck$key != "hgnc_symbol"))
  stop("[18] a ranked list is keyed by Ensembl id; every network target matches on HGNC symbol.")
emit(keycheck, "ranked_list_keycheck.csv")

# ---- the four named pre-2019 probes, resolved against the live vocabulary ----
probes <- bind_rows(lapply(TA$symbol_vocabulary_probes, function(p) {
  tibble(matrix_symbol  = p$matrix_symbol,
         current_symbol = p$current_symbol,
         matrix_symbol_in_ranked_treg  = p$matrix_symbol %in% names(STATS[[PRIMARY]]),
         current_symbol_in_ranked_treg = p$current_symbol %in% names(STATS[[PRIMARY]]),
         matrix_symbol_in_count_matrix  = p$matrix_symbol %in% MATRIX_SYMBOLS,
         current_symbol_in_count_matrix = p$current_symbol %in% MATRIX_SYMBOLS)
}))
print(as.data.frame(probes), row.names = FALSE)
emit(probes, "symbol_vocabulary_probes.csv")

# ============================================================================
# 3. The network and its variants
# ============================================================================

message("\n[3] Building the network variants ...")

net_raw <- readr::read_csv(NET_PATH, show_col_types = FALSE)
for (cc in c("source", "target", MOR_COL, SD_COL))
  if (!cc %in% names(net_raw))
    stop("[18] CollecTRI table lacks column ", cc)

net_base <- net_raw %>%
  transmute(source, target,
            mor = as.numeric(.data[[MOR_COL]]),
            sign_decision = .data[[SD_COL]]) %>%
  distinct(source, target, .keep_all = TRUE)
message(sprintf("  CollecTRI: %d unique edges, %d TFs, %d repressing edges",
                nrow(net_base), length(unique(net_base$source)), sum(net_base$mor < 0)))

# ---- HGNC alias resolution: current network symbol -> the symbol the matrix carries ----
# The machinery is shared across stages: helpers/symbol_alias.R holds
# one ownership guard and one rejection ledger for every consumer in the compartment, and
# this stage's two published ledger tables are the regression test that the lift changed
# nothing (02_analysis/scripts/00_symbol_alias_validate.R).
#
# What it decides. org.Hs.eg.db is the arbiter, and a target is recoverable when it is
# absent from the ranked list, its symbol resolves to exactly one Entrez id, and exactly
# one alias of that same Entrez id is present in the ranked list. The hazard the guard
# exists for: many retired symbols were reassigned as the official symbol of a DIFFERENT
# gene. PGF carries the alias PIGF, and PIGF now names a GPI-anchor biosynthesis gene;
# THPO carries TPO, and TPO now names thyroid peroxidase. Accepting either attaches one
# gene's expression to another gene's regulon edge, so a candidate that is the official
# symbol of any other Entrez id is rejected, counted and reported.
#
# The map the helper returns carries EVERY candidate with its resolution, so the two
# withholding classes the private copy returned silently — a reference symbol with two
# aliases in the universe, and one that is ambiguous in org.Hs.eg.db — are now rows rather
# than a shrug. Only `accepted` is ever applied to an edge.

make_variants <- function(universe) {
  signed <- net_base %>% filter(target %in% universe)
  unsigned <- signed %>% mutate(mor = 1)
  lit <- net_base %>% filter(target %in% universe, sign_decision != SD_DEF)
  missing_sym <- sort(setdiff(unique(net_base$target), universe))
  amap <- build_alias_map(missing_sym, universe, db = org.Hs.eg.db::org.Hs.eg.db)
  accepted <- amap %>% filter(resolution == "accepted")
  withheld <- amap %>% filter(resolution != "accepted")
  recovered <- net_base %>%
    filter(!target %in% universe) %>%
    rename(reference_symbol = target) %>%
    inner_join(accepted, by = "reference_symbol") %>%
    transmute(source, target = matrix_symbol, mor, sign_decision,
              recovered_from = reference_symbol)
  alias_net <- bind_rows(signed %>% mutate(recovered_from = NA_character_), recovered) %>%
    distinct(source, target, .keep_all = TRUE)
  list(signed = signed, unsigned = unsigned, literature_signed = lit,
       alias_recovered = alias_net %>% select(source, target, mor, sign_decision),
       .alias_map = amap, .accepted = accepted, .recovered_edges = recovered,
       .n_rejected_ambiguous = sum(
         amap$resolution == "rejected_symbol_belongs_to_another_gene"),
       .rejected = withheld)
}

NETS <- lapply(names(POPS), function(tag) make_variants(names(STATS[[tag]])))
names(NETS) <- names(POPS)

alias_recovery <- bind_rows(lapply(names(POPS), function(tag) {
  r <- NETS[[tag]]$.recovered_edges
  if (!nrow(r)) return(NULL)
  r %>% transmute(population = POPS[[tag]], tf = source,
                  reference_symbol = recovered_from, matrix_symbol = target, mor,
                  resolution = "accepted", focus_tf = source %in% FOCUS)
}))
# The withheld resolutions are published beside the accepted ones, at symbol level, so the
# guard is auditable and so a withholding never has to be inferred from an absence.
alias_rejected <- bind_rows(lapply(names(POPS), function(tag) {
  rj <- NETS[[tag]]$.rejected
  if (is.null(rj) || !nrow(rj)) return(NULL)
  net_base %>% filter(target %in% rj$reference_symbol) %>%
    transmute(population = POPS[[tag]], tf = source, reference_symbol = target,
              matrix_symbol = rj$matrix_symbol[match(target, rj$reference_symbol)], mor,
              resolution = rj$resolution[match(target, rj$reference_symbol)],
              focus_tf = source %in% FOCUS)
}))
alias_recovery <- bind_rows(alias_recovery, alias_rejected)
if (!nrow(alias_recovery))
  alias_recovery <- tibble(population = character(), tf = character(),
                           reference_symbol = character(), matrix_symbol = character(),
                           mor = numeric(), resolution = character(), focus_tf = logical())
emit(alias_recovery, "alias_recovery.csv")
acc <- alias_recovery %>% filter(resolution == "accepted")
message(sprintf("  alias recovery: %d edges accepted over %d TFs (%d on a focus TF); %d edges withheld, %d of them because the candidate is now another gene's official symbol",
                nrow(acc), length(unique(acc$tf)), sum(acc$focus_tf),
                sum(alias_recovery$resolution != "accepted"),
                sum(alias_recovery$resolution ==
                      "rejected_symbol_belongs_to_another_gene")))
message("  accepted resolutions on the focus TFs:")
print(as.data.frame(acc %>% filter(focus_tf, population == POPS[[PRIMARY]]) %>%
                      select(tf, reference_symbol, matrix_symbol)), row.names = FALSE)
message("  withheld resolutions (a sample of the distinct symbol pairs):")
print(as.data.frame(alias_recovery %>% filter(resolution != "accepted") %>%
                      distinct(reference_symbol, matrix_symbol, resolution) %>% head(15)),
      row.names = FALSE)

variant_summary <- bind_rows(lapply(names(POPS), function(tag) {
  bind_rows(lapply(VARIANTS, function(v) {
    n <- NETS[[tag]][[v]]
    sz <- n %>% count(source, name = "size")
    tibble(population = POPS[[tag]], variant = v,
           n_edges = nrow(n), n_tfs = nrow(sz),
           n_tfs_ge_minsize = sum(sz$size >= MINSZ),
           n_repressing_edges = sum(n$mor < 0),
           median_regulon_size = median(sz$size[sz$size >= MINSZ]),
           hif1a_size = sum(n$source == "HIF1A"),
           nfkb1_size = sum(n$source == "NFKB1"))
  }))
}))
print(as.data.frame(variant_summary %>% filter(population == POPS[[PRIMARY]])), row.names = FALSE)
emit(variant_summary, "network_variants.csv")

# ---- per-focus-TF accounting: the cause behind every unmatched target ----
vocab <- bind_rows(lapply(names(POPS), function(tag) {
  universe <- names(STATS[[tag]])
  amap <- NETS[[tag]]$.accepted
  bind_rows(lapply(FOCUS, function(tf) {
    t <- unique(net_base$target[net_base$source == tf])
    matched <- intersect(t, universe)
    unmatched <- setdiff(t, universe)
    expr_filtered <- intersect(unmatched, MATRIX_SYMBOLS)
    absent <- setdiff(unmatched, MATRIX_SYMBOLS)
    rec <- intersect(unmatched, amap$reference_symbol)
    tibble(population = POPS[[tag]], tf = tf,
           n_targets_in_network = length(t),
           n_matched = length(matched),
           n_unmatched = length(unmatched),
           n_expression_filtered = length(expr_filtered),
           n_absent_from_count_matrix = length(absent),
           n_alias_recoverable = length(rec),
           alias_recoverable_symbols = paste(sort(rec), collapse = "/"))
  }))
}))
print(as.data.frame(vocab %>% filter(population == POPS[[PRIMARY]]) %>%
                      select(-alias_recoverable_symbols)), row.names = FALSE)
emit(vocab, "symbol_vocabulary_check.csv")

# ============================================================================
# 4. The activity runs: ULM primary, MLM and consensus as forensics
# ============================================================================
# ULM regresses the contrast statistic of every gene on the mode of regulation the network
# assigns it and reports the slope's t, so the score is a statistic over the expression of
# the genes the network assigns to a factor, and it inherits that regulon's size,
# composition, target sharing and edge signs.

message("\n[4] decoupleR over ", length(POPS), " populations x ", length(VARIANTS),
        " variants x ", length(METHODS), " methods ...")

run_one <- function(mat, net) {
  # One decouple() call, so ULM, MLM and the consensus see identical network filtering and
  # a rank difference between them is the estimator alone.
  dec <- decouple(mat = mat, network = net,
                  .source = "source", .target = "target",
                  statistics = CONS_ST,
                  args = setNames(lapply(CONS_ST, function(s) list(.mor = "mor", minsize = MINSZ)),
                                  CONS_ST),
                  consensus_score = TRUE, minsize = MINSZ, show_toy_call = FALSE)
  dec %>% filter(statistic %in% METHODS)
}

ACT <- bind_rows(lapply(names(POPS), function(tag) {
  s   <- STATS[[tag]]
  mat <- matrix(s, ncol = 1, dimnames = list(names(s), sprintf("%s_SFvsPB", POPS[[tag]])))
  bind_rows(lapply(VARIANTS, function(v) {
    net <- NETS[[tag]][[v]] %>% select(source, target, mor)
    sz  <- net %>% count(source, name = "regulon_size")
    res <- run_one(mat, net)
    res %>%
      left_join(sz, by = "source") %>%
      group_by(statistic) %>%
      mutate(padj = p.adjust(p_value, "BH"),
             rank = rank(-score, ties.method = "first"),
             n_tfs_scored = n()) %>%
      ungroup() %>%
      transmute(population = POPS[[tag]], variant = v, method = statistic,
                tf = source, score, p_value, padj,
                rank = as.integer(rank), n_tfs_scored,
                pct_rank = 100 * rank / n_tfs_scored,
                regulon_size, direction = ifelse(score > 0, "Up", "Down"))
  }))
}))
message(sprintf("  %d rows over %d configurations",
                nrow(ACT), nrow(distinct(ACT, population, variant, method))))
emit(ACT %>% arrange(population, variant, method, rank), "tf_activity_all.csv")

# ============================================================================
# 5. The committed unsigned-regulon fgsea rows, the ranks being calibrated
# ============================================================================

message("\n[5] Reading the committed TF_Targets fgsea rows ...")
stopifnot(file.exists(FG_TABLE))
FG <- readr::read_csv(FG_TABLE, show_col_types = FALSE) %>%
  filter(database == FG_DB) %>%
  group_by(population) %>%
  arrange(desc(nes), .by_group = TRUE) %>%
  mutate(rank = as.integer(row_number()), n_tfs_scored = n(),
         pct_rank = 100 * rank / n_tfs_scored) %>%
  ungroup() %>%
  transmute(population, variant = "unsigned_geneset", method = "fgsea",
            tf = pathway_name, score = nes, p_value = pvalue, padj = padj_pooled,
            rank, n_tfs_scored, pct_rank, regulon_size = set_size,
            direction = ifelse(nes > 0, "Up", "Down"))
message(sprintf("  %d rows over %d populations; the %s family holds %d regulons",
                nrow(FG), n_distinct(FG$population), POPS[[PRIMARY]],
                sum(FG$population == POPS[[PRIMARY]])))

# The denominator a rank is read against. The sweep applies its size cap to the raw
# CollecTRI regulon, before intersecting with the ranked list, so a TF whose raw target
# list exceeds gsea_max_size sits outside that family. Those are the network's most
# promiscuous regulons, so the exclusion carries size with it and gets its own table.
raw_size <- net_base %>% distinct(source, target) %>% count(source, name = "raw_targets")
present_size <- NETS[[PRIMARY]]$signed %>% count(source, name = "targets_present")
CAPPED <- raw_size %>%
  left_join(present_size, by = "source") %>%
  filter(raw_targets > MAXSZ) %>%
  left_join(ACT %>% filter(population == POPS[[PRIMARY]], variant == "signed",
                           method == "ulm") %>%
              select(source = tf, ulm_score = score, ulm_rank = rank, ulm_padj = padj),
            by = "source") %>%
  mutate(population = POPS[[PRIMARY]],
         in_fgsea_family = source %in% FG$tf[FG$population == POPS[[PRIMARY]]],
         gsea_max_size = MAXSZ) %>%
  arrange(ulm_rank) %>%
  select(population, tf = source, raw_targets, targets_present, gsea_max_size,
         in_fgsea_family, ulm_score, ulm_rank, ulm_padj)
message(sprintf("  %d regulons sit above the sweep's raw-size cap of %d and are absent from that family:",
                nrow(CAPPED), MAXSZ))
print(as.data.frame(CAPPED %>% select(-population, -gsea_max_size)), row.names = FALSE)
emit(CAPPED, "fgsea_family_size_cap.csv")

# ============================================================================
# 6. Forensic 1: the rank cascade
# ============================================================================
# Four network variants crossed with three estimators, plus the committed fgsea rank. The
# mouse anchor's Hif1a moved #1 to #12 on a network swap and #12 to #142 from ULM to MLM,
# so those two axes are the ones this table travels along.

message("\n[6] The rank cascade ...")

CASCADE <- bind_rows(ACT, FG) %>%
  filter(tf %in% FOCUS) %>%
  mutate(configuration = paste(variant, method, sep = " / ")) %>%
  arrange(population, tf, variant, method) %>%
  select(population, tf, variant, method, configuration, score, p_value, padj,
         rank, n_tfs_scored, pct_rank, regulon_size, direction)
emit(CASCADE, "hif1a_rank_cascade.csv")

message("  HIF1A, primary population:")
print(as.data.frame(CASCADE %>% filter(tf == "HIF1A", population == POPS[[PRIMARY]]) %>%
                      select(configuration, score, rank, n_tfs_scored, padj, regulon_size)),
      row.names = FALSE)
message("  rank spread per focus TF (signed and unsigned_geneset variants, all methods):")
print(as.data.frame(CASCADE %>% filter(population == POPS[[PRIMARY]]) %>%
                      group_by(tf) %>%
                      summarise(best_rank = min(rank), worst_rank = max(rank),
                                spread = max(rank) - min(rank), .groups = "drop") %>%
                      arrange(best_rank)), row.names = FALSE)

# ============================================================================
# 7. Forensic 2: target decomposition against target promiscuity
# ============================================================================
# Which genes carry the score, each with its moderated t. Run for HIF1A and for NFKB1, so
# a reader has a comparator on the same axis.

message("\n[7] Target decomposition ...")

pop_tag  <- PRIMARY
universe <- names(STATS[[pop_tag]])
signed   <- NETS[[pop_tag]]$signed
de_prim  <- DE[[pop_tag]]

# How many CollecTRI regulons claim each target, over the universe the scores were computed
# on. A target in many regulons carries information about all of them jointly.
regulons_per_target <- signed %>% distinct(source, target) %>%
  count(target, name = "n_regulons")

DECOMP_TBL <- bind_rows(lapply(DECOMP, function(tf) {
  signed %>% filter(source == tf) %>%
    left_join(tibble(target = universe, stat = as.numeric(STATS[[pop_tag]])), by = "target") %>%
    left_join(regulons_per_target, by = "target") %>%
    left_join(de_prim %>% select(target = gene_symbol, avg_expr, log2FoldChange, padj_gene = padj),
              by = "target") %>%
    mutate(tf = tf,
           n_other_regulons = as.integer(n_regulons - 1L),
           contrib = sign(mor) * stat,
           selective = n_regulons <= SEL_MAX,
           promiscuity_band = cut(n_regulons, breaks = c(0, BANDS, Inf),
                                  labels = c(paste0("<=", BANDS), paste0(">", max(BANDS))),
                                  right = TRUE)) %>%
    select(tf, target, mor, sign_decision, stat, contrib, n_regulons, n_other_regulons,
           selective, promiscuity_band, avg_expr, log2FoldChange, padj_gene)
})) %>% arrange(tf, desc(contrib))
DECOMP_TBL$population <- POPS[[pop_tag]]
emit(DECOMP_TBL, "target_decomposition.csv")

decomp_summary <- DECOMP_TBL %>%
  group_by(tf) %>% mutate(total_contrib = sum(contrib)) %>%
  group_by(tf, promiscuity_band) %>%
  summarise(n_targets = n(), sum_contrib = sum(contrib),
            mean_contrib = mean(contrib), mean_n_regulons = mean(n_regulons),
            pct_of_total_contrib = 100 * sum(contrib) / first(total_contrib),
            .groups = "drop") %>%
  arrange(tf, promiscuity_band)
emit(decomp_summary, "target_decomposition_summary.csv")

sel_split <- DECOMP_TBL %>% group_by(tf) %>%
  summarise(n_targets = n(),
            n_selective = sum(selective),
            pct_targets_selective = 100 * mean(selective),
            median_n_regulons = median(n_regulons),
            mean_n_regulons = mean(n_regulons),
            total_contrib = sum(contrib),
            selective_contrib = sum(contrib[selective]),
            pct_contrib_selective = 100 * sum(contrib[selective]) / sum(contrib),
            .groups = "drop")
message("  selective (claimed by this regulon alone) versus promiscuous share of the score:")
print(as.data.frame(sel_split), row.names = FALSE)
message("  top ten contributing targets per decomposed TF:")
print(as.data.frame(DECOMP_TBL %>% group_by(tf) %>% slice_head(n = 10) %>% ungroup() %>%
                      select(tf, target, stat, contrib, n_regulons, mor)), row.names = FALSE)

# ============================================================================
# 8. Forensic 3: regulon-size calibration
# ============================================================================
# Activity against size across every TF tested, with a size-conditional residual taken
# over the real regulons and two random-regulon nulls. A systematic size dependence on
# this contrast bounds every large-regulon TF here.

message("\n[8] Regulon-size calibration ...")

ulm_prim <- ACT %>% filter(population == POPS[[pop_tag]], variant == "signed", method == "ulm")
fg_prim  <- FG  %>% filter(population == POPS[[pop_tag]])

CAL <- ulm_prim %>%
  select(tf, regulon_size, ulm_score = score, ulm_rank = rank, ulm_padj = padj) %>%
  left_join(fg_prim %>% select(tf, fgsea_nes = score, fgsea_rank = rank,
                               fgsea_set_size = regulon_size, fgsea_padj = padj),
            by = "tf") %>%
  mutate(population = POPS[[pop_tag]], focus_tf = tf %in% FOCUS,
         log_size = log10(regulon_size))

# Size-conditional expectation read off the real regulons: what activity a regulon of this
# size gets on this contrast. The residual is the part of a score size leaves unaccounted.
fit_resid <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  out <- rep(NA_real_, length(y)); fitv <- rep(NA_real_, length(y))
  if (sum(ok) > 10) {
    lo <- stats::loess(y[ok] ~ x[ok], span = 0.75, degree = 1)
    fitv[ok] <- as.numeric(stats::predict(lo))
    out[ok]  <- y[ok] - fitv[ok]
  }
  list(fit = fitv, resid = out)
}
r1 <- fit_resid(CAL$log_size, CAL$ulm_score)
CAL$ulm_score_size_expected <- r1$fit
CAL$ulm_score_size_residual <- r1$resid
r2 <- fit_resid(log10(CAL$fgsea_set_size), CAL$fgsea_nes)
CAL$fgsea_nes_size_expected <- r2$fit
CAL$fgsea_nes_size_residual <- r2$resid
# Rank on the residual as well, so one row carries activity rank beside size-adjusted rank.
rank_desc <- function(v) {
  out <- rep(NA_integer_, length(v)); ok <- is.finite(v)
  out[ok] <- as.integer(rank(-v[ok], ties.method = "first")); out
}
CAL$ulm_score_size_residual_rank   <- rank_desc(CAL$ulm_score_size_residual)
CAL$fgsea_nes_size_residual_rank   <- rank_desc(CAL$fgsea_nes_size_residual)
CAL$n_tfs_ulm_ranked   <- sum(is.finite(CAL$ulm_score_size_residual))
CAL$n_tfs_fgsea_ranked <- sum(is.finite(CAL$fgsea_nes_size_residual))
emit(CAL %>% arrange(desc(ulm_score)), "regulon_size_calibration.csv")
message("  activity rank against size-residual rank, focus TFs:")
print(as.data.frame(CAL %>% filter(focus_tf) %>%
                      select(tf, regulon_size, ulm_score, ulm_rank,
                             ulm_score_size_residual, ulm_score_size_residual_rank,
                             fgsea_nes, fgsea_rank, fgsea_nes_size_residual,
                             fgsea_nes_size_residual_rank) %>%
                      arrange(fgsea_rank)), row.names = FALSE, digits = 3)

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
# A label-permuted contrast is the reference. The size-versus-activity structure that
# survives shuffling the gene labels is arithmetic, and the rest belongs to the contrast's
# own broad shift, which a bigger regulon samples more thoroughly.
perm_stats <- STATS[[pop_tag]]
names(perm_stats) <- sample(names(perm_stats))
perm_mat <- matrix(perm_stats, ncol = 1,
                   dimnames = list(names(perm_stats), "permuted"))
perm_ulm <- run_ulm(mat = perm_mat, net = signed %>% select(source, target, mor),
                    .source = "source", .target = "target", .mor = "mor",
                    minsize = MINSZ) %>%
  left_join(signed %>% count(source, name = "regulon_size"), by = "source")

SPEAR <- bind_rows(lapply(names(POPS), function(tag) {
  a <- ACT %>% filter(population == POPS[[tag]], variant == "signed", method == "ulm")
  f <- FG  %>% filter(population == POPS[[tag]])
  tibble(population = POPS[[tag]],
         spearman_ulm_score_vs_size   = sp(a$regulon_size, a$score),
         spearman_ulm_abs_score_vs_size = sp(a$regulon_size, abs(a$score)),
         n_tfs_ulm = nrow(a),
         spearman_fgsea_nes_vs_size   = sp(f$regulon_size, f$score),
         n_tfs_fgsea = nrow(f))
})) %>%
  bind_rows(tibble(population = "Treg (gene labels permuted)",
                   spearman_ulm_score_vs_size = sp(perm_ulm$regulon_size, perm_ulm$score),
                   spearman_ulm_abs_score_vs_size = sp(perm_ulm$regulon_size, abs(perm_ulm$score)),
                   n_tfs_ulm = nrow(perm_ulm),
                   spearman_fgsea_nes_vs_size = NA_real_, n_tfs_fgsea = NA_integer_))
print(as.data.frame(SPEAR), row.names = FALSE)
emit(SPEAR, "regulon_size_spearman.csv")

# ---- random-regulon nulls, matched on size, then on size and expression ----
message("  size-matched random-regulon nulls, ", NDRAW, " draws per focus TF ...")

expr_lookup <- de_prim %>% select(gene_symbol, avg_expr) %>%
  filter(gene_symbol %in% universe) %>% distinct(gene_symbol, .keep_all = TRUE)
expr_vec <- setNames(expr_lookup$avg_expr, expr_lookup$gene_symbol)
expr_dec <- setNames(as.integer(cut(rank(expr_vec, ties.method = "first"),
                                    breaks = NDEC, labels = FALSE)), names(expr_vec))
by_decile <- split(names(expr_dec), expr_dec)

draw_matched <- function(obs_targets, match_expression) {
  if (!match_expression) return(sample(universe, length(obs_targets)))
  # Draw within the observed set's expression-decile composition, so the null keeps that
  # regulon's expression profile and varies gene identity alone.
  want <- table(expr_dec[obs_targets])
  out <- unlist(lapply(names(want), function(d) {
    pool <- by_decile[[d]]
    k <- min(as.integer(want[[d]]), length(pool))
    if (k > 0) sample(pool, k) else character(0)
  }), use.names = FALSE)
  short <- length(obs_targets) - length(out)
  if (short > 0) out <- c(out, sample(setdiff(universe, out), short))
  out
}

# One TF and one match mode per decoupleR call. Pooling every draw into a single network
# would build a genes-by-draws dense design matrix, and batching holds the estimator and
# the universe identical while the working set stays small.
prim_mat <- matrix(STATS[[pop_tag]], ncol = 1,
                   dimnames = list(names(STATS[[pop_tag]]), "primary"))

null_scores <- function(tf, match_expression) {
  obs <- signed %>% filter(source == tf)
  nn <- data.frame(
    source = rep(sprintf("d%05d", seq_len(NDRAW)), each = nrow(obs)),
    target = unlist(lapply(seq_len(NDRAW),
                           function(i) draw_matched(obs$target, match_expression)),
                    use.names = FALSE),
    # the observed sign composition is kept, so the null varies gene identity alone
    mor = unlist(lapply(seq_len(NDRAW), function(i) sample(obs$mor)), use.names = FALSE),
    stringsAsFactors = FALSE)
  u <- run_ulm(mat = prim_mat, net = nn, .source = "source", .target = "target",
               .mor = "mor", minsize = MINSZ)
  fg <- suppressWarnings(fgsea::fgsea(pathways = split(nn$target, nn$source),
                                      stats = STATS[[pop_tag]],
                                      minSize = MINSZ, maxSize = MAXSZ, nproc = 1))
  list(collectri_ulm_score = u$score, unsigned_geneset_fgsea_nes = fg$NES)
}

summarise_null <- function(vals, obs, tf, statistic, null_match) {
  vals <- vals[is.finite(vals)]
  tibble(tf = tf, statistic = statistic, null_match = null_match, obs = obs,
         n_draws = length(vals), null_mean = mean(vals), null_sd = stats::sd(vals),
         null_q95 = as.numeric(stats::quantile(vals, 0.95)), null_max = max(vals),
         pct_of_null = 100 * mean(vals < obs),
         n_null_ge_obs = sum(vals >= obs),
         p_empirical = (sum(vals >= obs) + 1) / (length(vals) + 1),
         z_vs_null = (obs - mean(vals)) / stats::sd(vals))
}

obs_u <- setNames(ulm_prim$score, ulm_prim$tf)
obs_f <- setNames(fg_prim$score,  fg_prim$tf)

NULLS <- bind_rows(lapply(c("size", "size_and_expression"), function(mode) {
  mx <- identical(mode, "size_and_expression")
  bind_rows(lapply(FOCUS, function(tf) {
    ns <- null_scores(tf, mx)
    message(sprintf("    %-6s %-20s ulm null mean %+0.2f q95 %+0.2f | obs %+0.2f",
                    tf, mode, mean(ns$collectri_ulm_score, na.rm = TRUE),
                    stats::quantile(ns$collectri_ulm_score, 0.95, na.rm = TRUE),
                    obs_u[[tf]] %||% NA_real_))
    bind_rows(
      summarise_null(ns$collectri_ulm_score, obs_u[[tf]] %||% NA_real_, tf,
                     "collectri_ulm_score", mode),
      if (!is.null(obs_f[[tf]]))
        summarise_null(ns$unsigned_geneset_fgsea_nes, obs_f[[tf]], tf,
                       "unsigned_geneset_fgsea_nes", mode)
    )
  }))
}))
NULLS <- NULLS %>%
  left_join(ulm_prim %>% select(tf, regulon_size), by = "tf") %>%
  mutate(population = POPS[[pop_tag]]) %>%
  select(population, tf, regulon_size, statistic, null_match, obs, n_draws,
         null_mean, null_sd, null_q95, null_max, pct_of_null, n_null_ge_obs,
         p_empirical, z_vs_null) %>%
  arrange(statistic, null_match, desc(obs))
print(as.data.frame(NULLS %>% filter(null_match == "size_and_expression")), row.names = FALSE)
emit(NULLS, "size_matched_null.csv")

# ============================================================================
# 9. Forensic 4: direction
# ============================================================================

message("\n[9] The direction audit ...")

# (a) What the recorded signs do. The committed TF_Targets rows pool activating and
#     repressing edges into one unsigned set, so those ranks use no sign. Compare signed
#     against unsigned within decoupleR, on one estimator.
sv <- ACT %>% filter(population == POPS[[pop_tag]], method == "ulm",
                     variant %in% c("signed", "unsigned")) %>%
  select(tf, variant, score, rank, padj, regulon_size) %>%
  pivot_wider(names_from = variant, values_from = c(score, rank, padj, regulon_size)) %>%
  left_join(signed %>% group_by(source) %>%
              summarise(n_repressing_edges = sum(mor < 0),
                        pct_repressing_edges = 100 * mean(mor < 0), .groups = "drop"),
            by = c("tf" = "source")) %>%
  mutate(population = POPS[[pop_tag]],
         delta_score = score_signed - score_unsigned,
         delta_rank  = rank_signed - rank_unsigned)
SIGNED_VS <- sv %>% filter(tf %in% FOCUS) %>%
  select(population, tf, regulon_size_signed, n_repressing_edges, pct_repressing_edges,
         score_signed, score_unsigned, delta_score,
         rank_signed, rank_unsigned, delta_rank, padj_signed, padj_unsigned) %>%
  arrange(rank_signed)
print(as.data.frame(SIGNED_VS), row.names = FALSE)
emit(SIGNED_VS, "signed_vs_unsigned.csv")
message(sprintf("  across all %d TFs: Spearman(signed rank, unsigned rank) = %.4f; %d of %d TFs move by >10 places",
                nrow(sv), sp(sv$rank_signed, sv$rank_unsigned),
                sum(abs(sv$delta_rank) > 10, na.rm = TRUE), nrow(sv)))

# (b) The canonical HIF1A-selective targets, by name, with their direction here. Every one
#     is reported with the cause when it fails to match, and with an alias probe, since a
#     silently dropped target is the failure mode under investigation.
alias_of <- function(sym) {
  m <- NETS[[pop_tag]]$.accepted
  hit <- m$matrix_symbol[m$reference_symbol == sym]
  if (length(hit)) hit[1] else NA_character_
}
CANON_TBL <- bind_rows(lapply(CANON, function(g) {
  in_rank  <- g %in% universe
  in_reg   <- g %in% signed$target[signed$source == "HIF1A"]
  mor_g    <- if (in_reg) signed$mor[signed$source == "HIF1A" & signed$target == g][1] else NA_real_
  de_row   <- de_prim %>% filter(gene_symbol == g)
  tibble(population = POPS[[pop_tag]], target = g,
         in_ranked_list = in_rank,
         in_count_matrix = g %in% MATRIX_SYMBOLS,
         in_hif1a_regulon = in_reg,
         mor = mor_g,
         stat = if (in_rank) as.numeric(STATS[[pop_tag]][[g]]) else NA_real_,
         log2FoldChange = if (nrow(de_row)) de_row$log2FoldChange[1] else NA_real_,
         padj_gene = if (nrow(de_row)) de_row$padj[1] else NA_real_,
         avg_expr = if (nrow(de_row)) de_row$avg_expr[1] else NA_real_,
         alias_in_matrix = if (in_rank) NA_character_ else alias_of(g),
         unmatched_cause = if (in_rank) NA_character_
                           else if (g %in% MATRIX_SYMBOLS) "expression_filter"
                           else "absent_from_count_matrix")
})) %>%
  mutate(direction_in_contrast = case_when(
    is.na(stat) ~ NA_character_, stat > 0 ~ "up_in_synovial_fluid",
    TRUE ~ "up_in_blood")) %>%
  arrange(desc(stat))
print(as.data.frame(CANON_TBL %>% select(target, in_ranked_list, in_hif1a_regulon, mor,
                                         stat, padj_gene, direction_in_contrast,
                                         unmatched_cause)), row.names = FALSE)
emit(CANON_TBL, "canonical_hif1a_targets.csv")

# ============================================================================
# 10. Checkpoints and the reproduction check
# ============================================================================

message("\n[10] Checkpoints and the reproduction check ...")
saveRDS(NETS, file.path(DIR_OBJ, "18_tf_networks.rds"))
saveRDS(list(activity = ACT, fgsea = FG, cascade = CASCADE, calibration = CAL,
             nulls = NULLS, decomposition = DECOMP_TBL, source_sha256 = NET_SHA),
        file.path(DIR_OBJ, "18_tf_activity.rds"))

# The committed sweep's TF_Targets rows are this stage's starting point, so the eight values
# it starts from are re-read here from the published table on every run.
gate <- FG %>% filter(population == POPS[[PRIMARY]], tf %in% FOCUS) %>%
  select(tf, nes = score, padj_pooled = padj, targets_tested = regulon_size, rank) %>%
  arrange(rank)
print(as.data.frame(gate), row.names = FALSE)

stopifnot(
  nrow(CASCADE) > 0,
  all(FOCUS %in% CASCADE$tf),
  # HIF1A's committed targets_tested equals the size its regulon has on this universe
  gate$targets_tested[gate$tf == "HIF1A"] ==
    sum(NETS[[PRIMARY]]$signed$source == "HIF1A"),
  nrow(DECOMP_TBL) > 0,
  all(is.finite(NULLS$obs))
)
message("\n[DONE] 18_tf_activity COMPUTE complete. Run 18_tf_activity_viz.R for figures.")
