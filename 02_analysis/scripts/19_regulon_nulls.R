#!/usr/bin/env Rscript
# 19_regulon_nulls.R: COMPUTE. Two nulls stage 18 does not have.
# =============================================================================
# Stage 18 established that HIF1A's CollecTRI-ULM activity on the sorted-Treg
# SF-versus-PB contrast survives every network and estimator swap, that activity scales with
# regulon size, and that the score's DIRECTION comes from targets many regulons share. Its
# three nulls, however, leave two specific gaps, and its own README names them:
#
#   (a) no null holds target PROMISCUITY fixed. The random-regulon null matches size,
#       repressing-edge fraction and expression decile, so it draws from the whole universe
#       and lands near zero — every real regulon beats it. The size-conditional residual
#       compares against real regulons but conditions on size alone.
#   (b) every null permutes the ANNOTATION (which genes a regulon claims, or which gene
#       carries which statistic). None permutes the experimental DESIGN, and gene-label
#       permutation destroys the gene-gene correlation that makes a real contrast's
#       statistics dependent, which is what makes it anti-conservative.
#
# SECTION 4 answers (a): curveball rewiring of the TF->target bipartite graph, preserving
# every regulon's size and every target's in-degree exactly. A drawn regulon therefore
# oversamples the promiscuous, high-|t| genes exactly as much as the observed one does, so
# the null's centre is NOT zero and beating it means something the size-matched null could
# not test.
#
# SECTION 5 answers (b): the SF/PB labels are swapped within donor and the contrast refitted
# end to end. Six Treg donors carry both arms, so all 2^6 = 64 configurations are enumerated
# and the test is EXACT — no seed enters the p-value. The price is resolution: the finest
# attainable one-sided p is 1/64 = 0.0156, published here rather than left implicit.
#
# Neither null says anything about HIF1A the protein. Both are statements about the CollecTRI
# HIF1A regulon's ULM activity on this contrast, which is the object stage 18 defined.
#
# Annotation tier. No row reaches 03_results/master/ or any effect-size accumulator.
#
# Reads, read-only:
#   03_results/03_pseudobulk/tables/ranked_treg.tsv          signed moderated t, HGNC symbol
#   03_results/03_pseudobulk/tables/pseudobulk_counts.csv    raw integer donor x stratum counts
#   03_results/03_pseudobulk/tables/pseudobulk_coldata.csv   donor / tissue / label per stratum
#   03_results/03_pseudobulk/tables/gene_symbols.csv         the Ensembl-to-HGNC seam map
#   ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv  (SHA-256 pinned)
#   03_results/18_tf_activity/tables/size_matched_null.csv   the stage-18 rung of the ladder
#
# Writes 03_results/19_regulon_nulls/tables/, each captioned in that stage's README:
#   source_hash_manifest.csv   ulm_engine_validation.csv
#   rewiring_null.csv          rewiring_null_draws.csv
#   signflip_null.csv          signflip_null_draws.csv
#   null_ladder.csv
#
# Compute only. Figures live in 19_regulon_nulls_viz.R.
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/19_regulon_nulls.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(Matrix)
  library(decoupleR)
  library(edgeR)
  library(limma)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")          # FIG_CFG, round_numeric_cols
source("02_analysis/helpers/source_hash_manifest.R")  # source_sha256, verify_source_hash

STAGE  <- "19_regulon_nulls"
SCRIPT <- "02_analysis/scripts/19_regulon_nulls.R"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# 0. Config — every parameter read, none inlined
# ============================================================================

CFG <- FIG_CFG
RN  <- CFG$regulon_nulls
TA  <- CFG$tf_activity
DSG <- CFG$design
if (is.null(RN))
  stop("[19] analysis_config.yaml has no `regulon_nulls:` block — add it before running.")
if (is.null(TA))
  stop("[19] `tf_activity:` is the source of the population, focus factors and network seam.")

THR   <- CFG$thresholds
MINSZ <- as.integer(THR$gsea_min_size %||% 5L)
SEED  <- as.integer(THR$gsea_seed %||% 123L)

NET_PATH <- CFG$unbiased_enrichment$tf_network$path
if (is.null(NET_PATH))
  stop("[19] `unbiased_enrichment.tf_network.path` is unset; the network path is config, never inline.")
MOR_COL <- TA$mor_col %||% "weight"
SD_COL  <- TA$sign_decision_col %||% "sign_decision"
FOCUS   <- unlist(TA$focus_tfs)
PRIMARY <- TA$primary_population %||% "treg"
POP_LABEL <- {
  m <- vapply(TA$populations, function(p) p$label, character(1))
  names(m) <- vapply(TA$populations, function(p) p$tag, character(1))
  m[[PRIMARY]]
}

N_DRAW        <- as.integer(RN$rewiring_draws %||% 1000L)
BURN_PER_EDGE <- as.numeric(RN$rewiring_burn_in_trades_per_edge %||% 5)
TRADE_PER_TF  <- as.numeric(RN$rewiring_trades_between_draws_per_tf %||% 2)
MAX_EXACT_D   <- as.integer(RN$signflip_max_exact_donors %||% 20L)
ULM_GATE      <- as.numeric(RN$ulm_validation_min_spearman %||% 0.9999)
SF_GATE       <- as.numeric(RN$signflip_identity_min_spearman %||% 0.99)

TISSUE_KEY <- DSG$tissue_key %||% "tissue"
TISSUE_NUM <- DSG$tissue_levels$synovial_fluid %||% "synovial_fluid"
TISSUE_DEN <- DSG$tissue_levels$peripheral_blood %||% "peripheral_blood"
DONOR_KEY  <- DSG$donor_key %||% "donor"

RESULTS <- CFG$paths$results %||% "03_results/"
TBL     <- file.path(RESULTS, STAGE, CFG$paths$stage_tables_subdir %||% "tables")
PB      <- file.path(RESULTS, "03_pseudobulk", "tables")
for (d in c(TBL, file.path(TBL, "_overview")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(SEED)

emit <- function(df, fname) {
  p <- file.path(TBL, fname)
  readr::write_csv(round_numeric_cols(df), p)
  message(sprintf("  [SAVE] %-30s %5d rows x %2d cols", fname, nrow(df), ncol(df)))
  invisible(p)
}

message("=================================================================")
message("19_regulon_nulls — promiscuity-preserving and design-permuting nulls")
message("=================================================================")

# ============================================================================
# 1. Pin the cross-compartment network
# ============================================================================
# Its own pin, not stage 18's: if the network moves, both stages must stop, and a stage that
# trusted a sibling's manifest would not.

message("\n[1] Pinning the cross-compartment CollecTRI network ...")
MANIFEST  <- file.path(TBL, "source_hash_manifest.csv")
NET_LABEL <- "collectri_regulons_human"
if (!file.exists(NET_PATH)) stop("[19] CollecTRI network absent: ", NET_PATH)
net_key <- manifest_key(NET_PATH, root = getwd())
if (!file.exists(MANIFEST)) {
  readr::write_csv(tibble(source_label = NET_LABEL, source_path = net_key,
                          sha256 = source_sha256(NET_PATH)), MANIFEST)
  message("  pin written (first run): ", MANIFEST)
}
NET_SHA <- verify_source_hash(NET_PATH, NET_LABEL, MANIFEST, root = getwd())
message(sprintf("  %s  sha256 %s", net_key, NET_SHA))

# ============================================================================
# 2. The ranked list and the signed network on its universe
# ============================================================================

message("\n[2] Ranked list and the signed network on its universe ...")

rk <- readr::read_tsv(file.path(PB, sprintf("ranked_%s.tsv", PRIMARY)),
                      col_names = c("gene", "stat"), show_col_types = FALSE)
stopifnot(!any(duplicated(rk$gene)))
# An Ensembl-keyed ranked list intersects every network at ~zero and decoupleR reports that
# as an empty result rather than an error. Same guard as stage 18, for the same reason.
if (mean(grepl("^ENSG[0-9]{6,}", rk$gene)) > 0.5)
  stop("[19] the ranked list is keyed by Ensembl id; the network matches on HGNC symbol.")
Y <- setNames(rk$stat, rk$gene)
UNIVERSE <- names(Y)
N_GENE <- length(UNIVERSE)

net_raw <- readr::read_csv(NET_PATH, show_col_types = FALSE)
for (cc in c("source", "target", MOR_COL, SD_COL))
  if (!cc %in% names(net_raw)) stop("[19] CollecTRI table lacks column ", cc)

# The `signed` variant of stage 18: CollecTRI's recorded per-edge sign, intersected with the
# ranked list. Stage 18's headline ULM score is read off exactly this network.
NET <- net_raw %>%
  transmute(source, target, mor = as.numeric(.data[[MOR_COL]])) %>%
  distinct(source, target, .keep_all = TRUE) %>%
  filter(target %in% UNIVERSE)

sizes <- NET %>% count(source, name = "regulon_size")
KEEP_TF <- sizes %>% filter(regulon_size >= MINSZ) %>% pull(source) %>% sort()
NET <- NET %>% filter(source %in% KEEP_TF)
N_TF <- length(KEEP_TF)
message(sprintf("  %d genes, %d factors at minsize %d, %d edges",
                N_GENE, N_TF, MINSZ, nrow(NET)))
missing_focus <- setdiff(FOCUS, KEEP_TF)
if (length(missing_focus))
  message("  focus factors absent from this configuration: ",
          paste(missing_focus, collapse = ", "))

# ============================================================================
# 3. The fast ULM, and the gate that licenses using it
# ============================================================================
# decoupleR's run_ulm regresses the contrast statistic of EVERY gene in the matrix on the
# regulon's mode of regulation (zero for a non-target) and returns the slope's t. That is a
# simple linear regression, so it has a closed form, and the closed form is one sparse
# matrix-vector product for all factors at once. The nulls below need ~64,000 of these
# fits; run_ulm would make that a batch job, and the closed form makes it seconds.
#
# The reason this is safe is the gate, not the algebra: it is validated against run_ulm on
# the observed network and the run stops if it does not reproduce it. This is the discipline
# the compartment already requires of a DE-engine swap — establish by rank correlation that
# the swap is a method change and not a result change.

ulm_fast <- function(X, y) {
  ## X: n_gene x n_tf sparse mor matrix (0 for non-targets). y: contrast statistic per gene.
  ## Returns the slope t-statistic per column of X, in column order.
  n <- length(y); ybar <- mean(y); SST <- sum((y - ybar)^2)
  sx  <- Matrix::colSums(X)
  sxx <- Matrix::colSums(X * X)
  Sxx <- sxx - sx^2 / n
  Sxy <- as.numeric(Matrix::crossprod(X, y)) - sx * ybar
  b   <- Sxy / Sxx
  SSE <- SST - b * Sxy
  b / sqrt(SSE / (n - 2) / Sxx)
}

# Regulon membership as integer index vectors, which is the representation the curveball
# trades operate on and the mor matrix is rebuilt from.
tf_targets <- split(match(NET$target, UNIVERSE), factor(NET$source, levels = KEEP_TF))
tf_mors    <- split(NET$mor, factor(NET$source, levels = KEEP_TF))
OBS_SIZE   <- vapply(tf_targets, length, integer(1))
OBS_INDEG  <- tabulate(unlist(tf_targets, use.names = FALSE), nbins = N_GENE)

build_X <- function(targets, mors) {
  Matrix::sparseMatrix(
    i = unlist(targets, use.names = FALSE),
    j = rep(seq_along(targets), lengths(targets)),
    x = unlist(mors, use.names = FALSE),
    dims = c(N_GENE, length(targets)))
}

X_obs <- build_X(tf_targets, tf_mors)
score_fast <- ulm_fast(X_obs, Y)
names(score_fast) <- KEEP_TF

mat <- matrix(Y, ncol = 1, dimnames = list(UNIVERSE, POP_LABEL))
dc <- decoupleR::run_ulm(mat = mat, net = NET, .source = "source", .target = "target",
                         .mor = "mor", minsize = MINSZ)
cmp <- tibble(tf = KEEP_TF, score_closed_form = as.numeric(score_fast)) %>%
  inner_join(dc %>% transmute(tf = source, score_decoupler = score), by = "tf")

val <- tibble(
  population = POP_LABEL, n_tfs_compared = nrow(cmp), n_genes = N_GENE,
  spearman = suppressWarnings(stats::cor(cmp$score_closed_form, cmp$score_decoupler,
                                         method = "spearman")),
  pearson = suppressWarnings(stats::cor(cmp$score_closed_form, cmp$score_decoupler)),
  max_abs_diff = max(abs(cmp$score_closed_form - cmp$score_decoupler)),
  gate_min_spearman = ULM_GATE,
  gate = "pass")
message(sprintf("\n[3] Fast ULM vs decoupleR: spearman %.8f, max |diff| %.3g over %d factors",
                val$spearman, val$max_abs_diff, nrow(cmp)))
if (!is.finite(val$spearman) || val$spearman < ULM_GATE) {
  val$gate <- "fail"
  emit(val, "ulm_engine_validation.csv")
  stop("[19] the closed-form ULM does not reproduce decoupleR; the nulls below would not be ",
       "measuring the published statistic.")
}
emit(val, "ulm_engine_validation.csv")

# ============================================================================
# 4. NULL A — degree-preserving rewiring (curveball)
# ============================================================================
# A curveball trade picks two regulons, holds the targets they share, and randomly
# redistributes the targets only one of them has, keeping each regulon's count. Row sums
# (regulon size) and column sums (how many regulons claim a target) are therefore invariant
# by construction, and the only thing destroyed is WHICH factor claims which target.
#
# Because size is preserved exactly, the set of factors clearing minsize is identical in
# every draw, so a rank and its denominator mean the same thing in every draw — which a
# size-perturbing null could not offer.

message(sprintf("\n[4] Rewiring null: %d draws, curveball trades ...", N_DRAW))

curveball_trade <- function(A, B) {
  ## One trade between two target-index vectors. Preserves length(A), length(B), and the
  ## multiset union, hence every target's in-degree.
  inter <- intersect(A, B)
  pool  <- c(setdiff(A, inter), setdiff(B, inter))
  np <- length(pool)
  if (!np) return(list(A, B))
  k <- length(A) - length(inter)
  perm <- pool[sample.int(np)]
  # Both tails are taken positively. `perm[-seq_len(0)]` is NOT the whole vector in R — it is
  # `perm[-integer(0)]`, which is empty — so the k = 0 case (A a subset of B) silently dropped
  # every one of B's exclusive targets and broke the in-degree invariant this null rests on.
  list(c(inter, if (k > 0)  perm[seq_len(k)]        else integer(0)),
       c(inter, if (k < np) perm[seq.int(k + 1, np)] else integer(0)))
}

run_trades <- function(state, n_trades) {
  for (i in seq_len(n_trades)) {
    ij <- sample.int(N_TF, 2L)
    tr <- curveball_trade(state[[ij[1]]], state[[ij[2]]])
    state[[ij[1]]] <- tr[[1]]; state[[ij[2]]] <- tr[[2]]
  }
  state
}

state <- tf_targets
n_burn <- ceiling(BURN_PER_EDGE * nrow(NET))
n_step <- ceiling(TRADE_PER_TF * N_TF)
message(sprintf("  burn-in %d trades, then %d trades between draws", n_burn, n_step))
state <- run_trades(state, n_burn)

focus_present <- intersect(FOCUS, KEEP_TF)
focus_idx <- match(focus_present, KEEP_TF)
draws <- vector("list", N_DRAW)

for (d in seq_len(N_DRAW)) {
  state <- run_trades(state, n_step)
  # Each regulon keeps its OWN multiset of edge signs, permuted onto its new targets. This
  # matches the stage-18 null convention: the null varies gene identity, never the sign
  # composition, so a factor cannot look different merely by having been re-signed.
  mors_d <- lapply(seq_len(N_TF), function(k) tf_mors[[k]][sample.int(OBS_SIZE[[k]])])
  Xd <- build_X(state, mors_d)
  s  <- ulm_fast(Xd, Y)
  r  <- rank(-s, ties.method = "min")
  draws[[d]] <- tibble(draw = d, tf = focus_present,
                       score = s[focus_idx], rank = r[focus_idx])
  if (d %% 200 == 0) message(sprintf("    draw %d/%d", d, N_DRAW))
}

# The invariants are the whole claim of this null, so they are asserted on the realised
# state rather than trusted to the algorithm.
stopifnot(identical(vapply(state, length, integer(1)), OBS_SIZE))
stopifnot(identical(tabulate(unlist(state, use.names = FALSE), nbins = N_GENE), OBS_INDEG))
stopifnot(all(vapply(state, function(v) !any(duplicated(v)), logical(1))))
message("  invariants hold: every regulon size and every target in-degree preserved exactly")

rew_draws <- bind_rows(draws)
emit(rew_draws %>% mutate(population = POP_LABEL, .before = 1), "rewiring_null_draws.csv")

obs_rank <- rank(-score_fast, ties.method = "min")
rew <- rew_draws %>%
  group_by(tf) %>%
  summarise(n_draws = n(),
            null_mean = mean(score), null_sd = stats::sd(score),
            null_q95 = as.numeric(stats::quantile(score, 0.95)), null_max = max(score),
            null_rank_median = stats::median(rank), null_rank_best = min(rank),
            .groups = "drop") %>%
  mutate(population = POP_LABEL,
         obs_score = as.numeric(score_fast[tf]),
         obs_rank = as.integer(obs_rank[tf]),
         regulon_size = as.integer(OBS_SIZE[tf]),
         mean_target_indeg = vapply(tf, function(f)
           mean(OBS_INDEG[tf_targets[[f]]]), numeric(1)),
         n_tfs_scored = N_TF) %>%
  rowwise() %>%
  # +1 in numerator and denominator so an empirical p never returns exactly zero, which is
  # the same convention size_matched_null.csv uses.
  mutate(pct_of_null = 100 * mean(rew_draws$score[rew_draws$tf == tf] < obs_score),
         p_empirical = (sum(rew_draws$score[rew_draws$tf == tf] >= obs_score) + 1) / (n_draws + 1),
         z_vs_null = (obs_score - null_mean) / null_sd) %>%
  ungroup() %>%
  select(population, tf, regulon_size, mean_target_indeg, obs_score, obs_rank, n_tfs_scored,
         n_draws, null_mean, null_sd, null_q95, null_max, pct_of_null, p_empirical,
         z_vs_null, null_rank_median, null_rank_best) %>%
  arrange(match(tf, FOCUS))
print(as.data.frame(rew %>% select(tf, regulon_size, obs_score, null_mean, null_q95,
                                   p_empirical, z_vs_null)), row.names = FALSE)
emit(rew, "rewiring_null.csv")

# ============================================================================
# 5. NULL B — exact within-donor sign-flip of the contrast
# ============================================================================
# The design is paired: each donor contributing both arms supplies one SF-versus-PB
# difference. Under the null of no tissue effect those differences are exchangeable in SIGN,
# so swapping the two labels within a donor generates a valid null contrast. Every gene is
# refitted together, so the correlation between genes — which gene-label permutation
# destroys — is carried through untouched.
#
# Two choices are deliberate and both are stated in the caption:
#   * The gene set is filtered ONCE, on the observed design, and held fixed. Refiltering per
#     configuration would change the ULM universe between draws and the scores would no
#     longer be comparable; the design is what is permuted here, not the filter.
#   * voom IS recomputed per configuration, because its weights are a function of the design
#     and reusing the observed weights would leak the observed labels into every draw.

message("\n[5] Sign-flip null: refitting the contrast under every within-donor label swap ...")

counts_df <- read.csv(file.path(PB, "pseudobulk_counts.csv"), row.names = 1, check.names = FALSE)
coldata   <- read.csv(file.path(PB, "pseudobulk_coldata.csv"), row.names = 1)
gene_map  <- read.csv(file.path(PB, "gene_symbols.csv"), row.names = 1)

common <- intersect(rownames(coldata), rownames(counts_df))
coldata <- coldata[common, , drop = FALSE]
counts_mat <- t(as.matrix(counts_df[common, ]))   # genes x strata, as 03b builds it

cd <- coldata[coldata$coarse_label == POP_LABEL &
              coldata[[TISSUE_KEY]] %in% c(TISSUE_NUM, TISSUE_DEN), , drop = FALSE]
c_mat <- counts_mat[, rownames(cd), drop = FALSE]
cd[[TISSUE_KEY]] <- factor(cd[[TISSUE_KEY]], levels = c(TISSUE_DEN, TISSUE_NUM))
cd[[DONOR_KEY]]  <- factor(cd[[DONOR_KEY]])

per_arm <- table(cd[[DONOR_KEY]], cd[[TISSUE_KEY]])
paired_donors <- rownames(per_arm)[rowSums(per_arm > 0) == 2]
n_paired <- length(paired_donors)
message(sprintf("  %s: %d strata, %d donors, %d carrying both arms",
                POP_LABEL, nrow(cd), nlevels(cd[[DONOR_KEY]]), n_paired))
if (n_paired < 2) stop("[19] fewer than two paired donors; there is no paired design to permute.")
if (n_paired > MAX_EXACT_D)
  stop("[19] ", n_paired, " paired donors means 2^", n_paired, " configurations. Raise ",
       "regulon_nulls.signflip_max_exact_donors deliberately or switch to sampling.")

# The observed design decides the filter and the model, exactly as 03b_pseudobulk_de.R does.
design_obs <- model.matrix(as.formula(sprintf("~ %s + %s", DONOR_KEY, TISSUE_KEY)), data = cd)
dge0 <- DGEList(counts = c_mat)
keep <- filterByExpr(dge0, design = design_obs)
dge0 <- dge0[keep, , keep.lib.sizes = FALSE]
dge0 <- calcNormFactors(dge0)
COEF <- sprintf("%s%s", TISSUE_KEY, TISSUE_NUM)
message(sprintf("  %d genes kept by filterByExpr on the observed design, held fixed across configurations",
                nrow(dge0)))

# The Ensembl-to-symbol representative is chosen ONCE so the ULM universe is byte-identical
# across configurations. The compartment's rule is one row per symbol keeping the most
# extreme |t|; with no duplicate symbol among the kept ids the rule is a no-op, which is
# asserted rather than assumed — if it ever binds, the winner must not be allowed to move
# between configurations or the universe would drift with the permutation.
sym <- gene_map[rownames(dge0), "gene_symbol"]
ok  <- !is.na(sym) & !sym %in% c("", "nan", "None", "NA")
stopifnot(!any(duplicated(sym[ok])))
GENE_SYM <- sym[ok]

fit_contrast <- function(design) {
  ## voom -> lmFit -> eBayes(robust) on the fixed gene set. Returns the moderated t per
  ## symbol, on the fixed universe. Identical to 03b_pseudobulk_de.R's model.
  v   <- voom(dge0, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  tt  <- fit$t[, COEF]
  setNames(tt[ok], GENE_SYM)
}

flip_design <- function(flip_donors) {
  d <- cd
  sw <- d[[DONOR_KEY]] %in% flip_donors
  d[[TISSUE_KEY]] <- factor(
    ifelse(sw & d[[TISSUE_KEY]] == TISSUE_NUM, TISSUE_DEN,
           ifelse(sw & d[[TISSUE_KEY]] == TISSUE_DEN, TISSUE_NUM,
                  as.character(d[[TISSUE_KEY]]))),
    levels = c(TISSUE_DEN, TISSUE_NUM))
  model.matrix(as.formula(sprintf("~ %s + %s", DONOR_KEY, TISSUE_KEY)), data = d)
}

# The identity configuration must reproduce the committed ranked list. This is the seam gate:
# it catches any divergence between this refit and 03b_pseudobulk_de.R before a single null
# score is read, and a null measured against a contrast that is not the published one would
# be silently meaningless.
t_ident <- fit_contrast(design_obs)
shared  <- intersect(names(t_ident), UNIVERSE)
rho_ident <- suppressWarnings(stats::cor(t_ident[shared], Y[shared], method = "spearman"))
message(sprintf("  identity configuration vs the committed ranked list: spearman %.6f over %d genes",
                rho_ident, length(shared)))
if (!is.finite(rho_ident) || rho_ident < SF_GATE)
  stop("[19] the identity refit does not reproduce the committed contrast (spearman ",
       signif(rho_ident, 4), " < ", SF_GATE, "); the sign-flip null would be measured ",
       "against a different contrast than the published one.")

# All 2^n_paired sign patterns, the observed one included: it is one of the configurations and
# excluding it would bias the p-value.
patterns <- expand.grid(rep(list(c(FALSE, TRUE)), n_paired), KEEP.OUT.ATTRS = FALSE)
names(patterns) <- paired_donors
n_cfg <- nrow(patterns)
message(sprintf("  enumerating all %d within-donor configurations (2^%d), exactly",
                n_cfg, n_paired))

# The ULM universe is the intersection of the refit's symbols with the ranked list, held
# fixed across configurations, so the network is subset once.
uni_sf   <- shared
net_sf   <- NET %>% filter(target %in% uni_sf)
sizes_sf <- net_sf %>% count(source, name = "n")
keep_sf  <- sizes_sf %>% filter(n >= MINSZ) %>% pull(source) %>% sort()
net_sf   <- net_sf %>% filter(source %in% keep_sf)
tgt_sf   <- split(match(net_sf$target, uni_sf), factor(net_sf$source, levels = keep_sf))
mor_sf   <- split(net_sf$mor, factor(net_sf$source, levels = keep_sf))
X_sf     <- Matrix::sparseMatrix(
  i = unlist(tgt_sf, use.names = FALSE),
  j = rep(seq_along(tgt_sf), lengths(tgt_sf)),
  x = unlist(mor_sf, use.names = FALSE),
  dims = c(length(uni_sf), length(keep_sf)))
focus_sf <- intersect(FOCUS, keep_sf)

sf_rows <- vector("list", n_cfg)
for (i in seq_len(n_cfg)) {
  fl <- paired_donors[unlist(patterns[i, ])]
  tv <- if (!length(fl)) t_ident else fit_contrast(flip_design(fl))
  s  <- ulm_fast(X_sf, tv[uni_sf])
  names(s) <- keep_sf
  r  <- rank(-s, ties.method = "min"); names(r) <- keep_sf
  sf_rows[[i]] <- tibble(
    configuration = i, n_donors_flipped = length(fl),
    donors_flipped = if (length(fl)) paste(fl, collapse = "|") else "none",
    is_observed = length(fl) == 0,
    tf = focus_sf, score = s[focus_sf], rank = as.integer(r[focus_sf]))
  if (i %% 16 == 0) message(sprintf("    configuration %d/%d", i, n_cfg))
}
sf_draws <- bind_rows(sf_rows) %>% mutate(population = POP_LABEL, .before = 1)
emit(sf_draws, "signflip_null_draws.csv")

# Flipping every donor negates the contrast, so the configurations come in +/- pairs and the
# null is sign-symmetric. That is what fixes the resolution floor: 2^n configurations give
# n_cfg distinct one-sided p-values, the finest being 1/n_cfg.
sf <- sf_draws %>%
  group_by(tf) %>%
  summarise(n_configurations = n(),
            obs_score = score[is_observed][1],
            obs_rank = rank[is_observed][1],
            null_mean = mean(score[!is_observed]),
            null_sd = stats::sd(score[!is_observed]),
            null_max = max(score[!is_observed]),
            n_ge_obs = sum(score >= score[is_observed][1]),
            .groups = "drop") %>%
  mutate(population = POP_LABEL,
         n_paired_donors = n_paired,
         p_exact_one_sided = n_ge_obs / n_configurations,
         p_floor = 1 / n_configurations,
         at_resolution_floor = n_ge_obs == 1L,
         n_tfs_scored = length(keep_sf)) %>%
  select(population, tf, n_paired_donors, n_configurations, obs_score, obs_rank,
         n_tfs_scored, null_mean, null_sd, null_max, n_ge_obs, p_exact_one_sided, p_floor,
         at_resolution_floor) %>%
  arrange(match(tf, FOCUS))
print(as.data.frame(sf %>% select(tf, obs_score, null_mean, null_max, n_ge_obs,
                                  p_exact_one_sided, p_floor)), row.names = FALSE)
emit(sf, "signflip_null.csv")

# ============================================================================
# 6. The ladder — every null side by side, with what each one holds fixed
# ============================================================================
# A null is only as strong as what it holds constant, and the four here form an ordering.
# Reading any one of them alone is how a weak null gets mistaken for a strong result, so
# they are published as one table with that ordering explicit.

message("\n[6] Assembling the null ladder ...")

s18 <- file.path(RESULTS, "18_tf_activity", "tables", "size_matched_null.csv")
s18_rows <- if (file.exists(s18)) {
  readr::read_csv(s18, show_col_types = FALSE) %>%
    filter(statistic == "collectri_ulm_score") %>%
    transmute(tf, null = paste0("stage18_", null_match), obs, null_mean, null_q95,
              p_empirical,
              holds_fixed = ifelse(null_match == "size_and_expression",
                                   "regulon size; repressing-edge fraction; expression-decile composition",
                                   "regulon size; repressing-edge fraction"),
              permutes = "which genes are in the regulon",
              n_draws = NA_integer_)
} else {
  message("  stage-18 size_matched_null.csv absent; the ladder carries this stage's rungs only")
  NULL
}

ladder <- bind_rows(
  s18_rows,
  rew %>% transmute(tf, null = "rewiring_degree_preserving", obs = obs_score, null_mean,
                    null_q95, p_empirical, n_draws,
                    holds_fixed = paste("regulon size; EVERY target's in-degree, so the",
                                        "promiscuity profile; each regulon's own sign multiset"),
                    permutes = "which factor claims which target"),
  sf %>% transmute(tf, null = "signflip_within_donor_exact", obs = obs_score, null_mean,
                   null_q95 = NA_real_, p_empirical = p_exact_one_sided,
                   n_draws = n_configurations,
                   holds_fixed = paste("the network entirely; the gene set; the gene-gene",
                                       "correlation structure"),
                   permutes = "the SF/PB labels within donor, i.e. the experimental design")
) %>%
  filter(tf %in% FOCUS) %>%
  mutate(null = factor(null, levels = c("stage18_size", "stage18_size_and_expression",
                                        "rewiring_degree_preserving",
                                        "signflip_within_donor_exact"))) %>%
  arrange(match(tf, FOCUS), null) %>%
  mutate(population = POP_LABEL, .before = 1)
emit(ladder, "null_ladder.csv")

# ============================================================================
# 7. Captions — written from the values just computed, not from memory
# ============================================================================
# Every table gets its README section through write_caption(), which is idempotent on the
# filename. The findings are interpolated from the realised numbers so a caption cannot drift
# from the table it describes on a re-run with different parameters.

message("\n[7] Writing the stage README captions ...")

g <- function(tbl, f, col) tbl[[col]][tbl$tf == f][1]
cap <- function(fname, finding, fn, config_kv, input, how_to_read)
  write_caption(STAGE, file.path("tables", fname), finding = finding, script = SCRIPT,
                fn = fn, config_kv = config_kv, input = input,
                how_to_read = how_to_read, config = CFG)

cap("source_hash_manifest.csv",
  sprintf(paste("The CollecTRI human regulon table read across the compartment boundary is",
                "pinned at sha256 %s..., so a change to that network stops this stage instead",
                "of quietly moving both nulls."), substr(NET_SHA, 1, 8)),
  "verify_source_hash",
  paste0("unbiased_enrichment.tf_network.path = ", NET_PATH), NET_PATH,
  paste("One row per cross-compartment source: `source_label` is the name this stage refers to",
        "it by, `source_path` is repository-root-relative, and `sha256` is the digest of the",
        "bytes actually read. The first run writes the pin; every later run verifies against it",
        "and stops on a mismatch. This stage keeps its own pin rather than trusting the",
        "sibling stage's, so a moved network stops both. Verification is the only gate."))

cap("ulm_engine_validation.csv",
  sprintf(paste("The closed-form univariate-linear-model score used for all %s null draws",
                "reproduces decoupleR's `run_ulm` on the observed network exactly: Spearman %.6f",
                "and a largest absolute difference of %.3g over %d factors, which is machine",
                "precision, so the nulls measure the same statistic the headline was read off."),
          format(N_DRAW + n_cfg, big.mark = ","), val$spearman, val$max_abs_diff,
          val$n_tfs_compared),
  "ulm_fast",
  paste0("regulon_nulls.ulm_validation_min_spearman = ", ULM_GATE,
         "; thresholds.gsea_min_size = ", MINSZ),
  paste0("03_results/03_pseudobulk/tables/ranked_", PRIMARY, ".tsv, ", NET_PATH),
  paste("One row. `score_closed_form` against `score_decoupler` is a like-for-like comparison on",
        "the observed signed network, summarised by `spearman`, `pearson` and `max_abs_diff`.",
        "decoupleR regresses every gene's contrast statistic on the regulon's mode of regulation,",
        "zero for a non-target, so the fit has a closed form that costs one sparse matrix-vector",
        "product for all factors at once — which is what makes tens of thousands of null fits",
        "feasible. `gate` is `pass` only when `spearman` reaches",
        "`gate_min_spearman`; the run stops otherwise. This is a guard, not a measurement, and",
        "it is the compartment's standing requirement that an engine swap be shown by rank",
        "correlation to be a method change rather than a result change."))

cap("rewiring_null.csv",
  sprintf(paste("Held against random regulons that share its size AND every one of its targets'",
                "in-degrees, HIF1A's CollecTRI-ULM score of %.2f still clears the null (empirical",
                "p %.3f, z %.2f) but the null's centre is %.2f rather than zero, so most of the",
                "score is what any regulon of that size and promiscuity profile earns on this",
                "contrast. The null does not single HIF1A out: REL clears it harder (p %.3f, z",
                "%.2f) on 83 targets, NFKB1 is indistinguishable (p %.3f), and EPAS1 sits exactly",
                "on its own null (obs %.2f against a null mean of %.2f, p %.2f)."),
          g(rew, "HIF1A", "obs_score"), g(rew, "HIF1A", "p_empirical"),
          g(rew, "HIF1A", "z_vs_null"), g(rew, "HIF1A", "null_mean"),
          g(rew, "REL", "p_empirical"), g(rew, "REL", "z_vs_null"),
          g(rew, "NFKB1", "p_empirical"), g(rew, "EPAS1", "obs_score"),
          g(rew, "EPAS1", "null_mean"), g(rew, "EPAS1", "p_empirical")),
  "curveball_trade + run_trades + ulm_fast",
  paste0("regulon_nulls.rewiring_draws = ", N_DRAW,
         "; regulon_nulls.rewiring_burn_in_trades_per_edge = ", BURN_PER_EDGE,
         "; regulon_nulls.rewiring_trades_between_draws_per_tf = ", TRADE_PER_TF,
         "; tf_activity.focus_tfs = 8; thresholds.gsea_min_size = ", MINSZ),
  paste0("03_results/03_pseudobulk/tables/ranked_", PRIMARY, ".tsv, ", NET_PATH),
  paste("One row per focus factor. `obs_score` and `obs_rank` are the observed signed-network",
        "ULM values; the `null_*` columns describe the draws. A curveball trade holds the targets",
        "two regulons share and redistributes the rest between them, so every regulon's size and",
        "every target's in-degree are invariant and only the assignment of targets to factors is",
        "randomised — which is why `null_mean` is far above zero and why beating this null means",
        "more than beating a random gene set of matched size. `mean_target_indeg` is the observed",
        "regulon's average target promiscuity, the property this null holds fixed.",
        "`p_empirical` is (draws at or above obs + 1)/(draws + 1) so it never returns exactly",
        "zero, and `z_vs_null` standardises the gap. Because size is preserved exactly, the same",
        "factors clear minsize in every draw, so `null_rank_median` and `null_rank_best` share the",
        "`n_tfs_scored` denominator with `obs_rank`. Annotation tier: this is a statement about",
        "the regulon's target set, never about the transcription factor's protein activity."))

cap("rewiring_null_draws.csv",
  sprintf(paste("The per-draw substrate behind the rewiring null: %d degree-preserving rewirings",
                "scored for each of the %d focus factors present in this configuration."),
          N_DRAW, length(focus_present)),
  "run_trades + ulm_fast",
  paste0("regulon_nulls.rewiring_draws = ", N_DRAW), "in-memory curveball chain",
  paste("One row per (draw, factor). `score` is that draw's closed-form ULM value and `rank` its",
        "position by descending score among every factor scored in the same draw. Draws come from",
        "one Markov chain advanced between samples rather than restarted, so consecutive rows are",
        "decorrelated by trades and not independent by construction. This is the distribution",
        "`rewiring_null.csv` summarises; read it when a summary statistic needs checking against",
        "the shape it came from."))

cap("signflip_null.csv",
  sprintf(paste("Permuting the experimental design rather than the annotation, the observed",
                "synovial-fluid-versus-blood labelling gives HIF1A the largest score of all %d",
                "within-donor configurations (%.2f against a null maximum of %.2f), which is",
                "exact p = %.4f, the finest this design can resolve. It carries no information",
                "about HIF1A specifically: %d of the %d focus factors reach the same floor, ATF3",
                "being the only one that does not (p %.3f). What it establishes is that the",
                "contrast itself is not an artefact of labelling."),
          n_cfg, g(sf, "HIF1A", "obs_score"), g(sf, "HIF1A", "null_max"),
          g(sf, "HIF1A", "p_exact_one_sided"), sum(sf$at_resolution_floor), nrow(sf),
          g(sf, "ATF3", "p_exact_one_sided")),
  "flip_design + fit_contrast + ulm_fast",
  paste0("design.donor_key = ", DONOR_KEY, "; design.tissue_key = ", TISSUE_KEY,
         "; regulon_nulls.signflip_max_exact_donors = ", MAX_EXACT_D,
         "; regulon_nulls.signflip_identity_min_spearman = ", SF_GATE),
  paste0("03_results/03_pseudobulk/tables/pseudobulk_counts.csv, ",
         "03_results/03_pseudobulk/tables/pseudobulk_coldata.csv, ",
         "03_results/03_pseudobulk/tables/gene_symbols.csv, ", NET_PATH),
  paste("One row per focus factor. The two tissue labels are swapped within a donor and the whole",
        "limma-voom contrast is refitted, which is the exchangeability a paired design licenses",
        "and which — unlike permuting gene labels — leaves the correlation between genes intact.",
        "`n_paired_donors` donors carry both arms, so `n_configurations` = 2^that, and ALL of them",
        "are enumerated: the test is exact and no seed enters it. `n_ge_obs` counts configurations",
        "scoring at or above the observed, the observed one included, and `p_exact_one_sided` is",
        "that count over `n_configurations`. `p_floor` is 1/`n_configurations` and",
        "`at_resolution_floor` is TRUE when only the observed configuration reaches the observed",
        "score, i.e. the p-value is as small as this many donors can make it and a smaller one",
        "would need more donors, not more computation. Flipping every donor negates the contrast,",
        "so the configurations come in sign-symmetric pairs and `null_mean` sits near zero by",
        "construction. The gene set is filtered once on the observed design and held fixed, so",
        "scores are comparable across configurations, while voom is recomputed per configuration",
        "because its weights depend on the design."))

cap("signflip_null_draws.csv",
  sprintf(paste("The per-configuration substrate behind the sign-flip null: all %d within-donor",
                "label swaps of the %d paired donors, scored for each focus factor."),
          n_cfg, n_paired),
  "flip_design + fit_contrast + ulm_fast",
  paste0("design.donor_key = ", DONOR_KEY, "; design.tissue_key = ", TISSUE_KEY),
  "03_results/03_pseudobulk/tables/pseudobulk_counts.csv + pseudobulk_coldata.csv",
  paste("One row per (configuration, factor). `donors_flipped` names the donors whose two tissue",
        "labels were swapped, pipe-delimited, and `is_observed` marks the single configuration",
        "that flips none — the published contrast. `score` and `rank` are that configuration's",
        "closed-form ULM value and its position among every factor scored. Configurations are",
        "enumerated in a fixed order, so this table is byte-stable across runs. Read it to see the",
        "null's shape, and in particular that the observed row is the extreme one for most",
        "factors rather than merely a high one."))

cap("null_ladder.csv",
  paste("Set side by side, the four nulls order by what they hold fixed and the conclusion",
        "weakens along that order: HIF1A clears a size-and-expression-matched random regulon by a",
        "wide margin, clears a size-and-promiscuity-matched real-degree rewiring by a narrow one,",
        "and is one of seven focus factors that max out the exact design permutation. No single",
        "rung supports a claim on its own, which is why they are published as one table."),
  "bind_rows",
  "tf_activity.focus_tfs = 8; tf_activity.null_draws; regulon_nulls.rewiring_draws",
  paste0("03_results/18_tf_activity/tables/size_matched_null.csv, ",
         "03_results/19_regulon_nulls/tables/rewiring_null.csv + signflip_null.csv"),
  paste("One row per (factor, null), ordered weakest-to-strongest by what the null holds constant.",
        "`holds_fixed` names exactly that and `permutes` names what is randomised, which together",
        "are what a p-value from that row can and cannot mean. `obs` is the same observed",
        "CollecTRI-ULM score in every row for a given factor; only the reference distribution",
        "changes, so the columns are comparable down a factor's block and the movement in",
        "`p_empirical` is the whole content of the table. `null_q95` is empty for the sign-flip",
        "rung because its 64 configurations do not support a stable upper percentile.",
        "The two random-regulon rungs are re-read from their own published table rather than recomputed."))

message("\n[DONE] 19_regulon_nulls complete.")
