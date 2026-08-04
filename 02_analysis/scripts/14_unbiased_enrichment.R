#!/usr/bin/env Rscript
# 14_unbiased_enrichment.R — COMPUTE (no plotting)
# =============================================================================
# What this stage asks, and why it is separate from stage 05
# -----------------------------------------------------------------------------
# Stage 05 asks a TARGETED question of the JIA niche contrast: does one named,
# mouse-derived signature enrich in synovial fluid relative to paired blood? A
# targeted test can only answer about the set it was handed. This stage asks the
# unbiased counterpart off the SAME frozen ranked lists: what does the
# synovial-fluid-versus-paired-blood contrast contain at all, across curated
# databases, with no set privileged? Only the two together tell a reader whether
# the mouse-derived enrichment is a distinctive feature of this contrast or one of
# many co-enriching inflammatory programs.
#
# Two methods, deliberately different in what they need:
#   (1) Pre-ranked fgsea over every set in seven human MSigDB collections, the
#       CollecTRI TF regulons, the toolkit's human MitoPathways build, this
#       compartment's frozen curated lists (the Hallmark re-pins, the curated
#       heat-shock-response lens and the curated TCR activation lens), the
#       mouse-derived projected UP arms, and the frozen SAVI axes. Benjamini-
#       Hochberg is applied BOTH per database (comparable to a single-collection
#       run) and POOLED across the whole family of tests within one population, so
#       a headline can be read against how many hypotheses were asked of the same
#       ranked list.
#   (2) decoupleR MLM on the human PROGENy model — fourteen signalling footprints
#       with continuous weights, needing no gene-set list at all, so it is not
#       exposed to the size and curation choices method (1) inherits.
#
# THE POOLED FAMILY CONTAINS EACH SET ONCE. Two collections in this sweep legitimately
# carry the same gene set: `project_frozen` pins six MSigDB Hallmark sets to files so
# the decomposition and purge stages have a size-validated asset that cannot move under
# a msigdbr upgrade, and the `Hallmark` collection then fetches those same six live.
# Scoring both is fine; letting both into one pooled Benjamini-Hochberg family is not,
# because the family then contains exact duplicate hypotheses and every population's
# rank denominator is six too high. Section 2b resolves the collision structurally, so
# a set id present in more than one collection is scored in both (each per-database
# table still reads as a standalone single-collection run) but pooled exactly once.
#
# CLAIM TIER. Every number here is a ranked-list enrichment statistic or a
# footprint activity score on the same donor-pseudobulk contrast the confirmatory
# spine already published. Nothing here creates a new claim, nothing is written to
# 03_results/master/, and language stays correlative: a set enriching says its gene
# content moves with the synovial-fluid side of this contrast, not that the program
# it is named for is present or is driving anything.
#
# THE CORRECTNESS GATE. WT_heat_up is run as an ordinary member of the sweep, and
# its NES must land on the published stage-05 value off the same ranked list. If it
# does not, the set prep, ranking or fgsea parameters differ from the published
# pipeline and NOTHING in the sweep is comparable to it — so the gate STOPS the run
# before the expensive collections are touched. It is never retuned to agree.
#
# THE SILENT-FAILURE MODE THIS GUARDS AGAINST. Count matrices in this compartment
# are keyed by Ensembl id; every reference gene set matches on HGNC symbol. A ranked
# list that leaked Ensembl ids intersects every collection at approximately zero,
# and fgsea reports that as empty/NA rather than as an error — it looks like a
# biological null. So the ranked lists are key-checked before anything runs, and the
# per-database overlap counts are published as a first-class table rather than
# assumed.
#
# Inputs (READ-ONLY):
#   03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv   signed moderated t, HGNC symbol
#   03_results/03_pseudobulk/tables/de_SFvsPB_{treg,tcon,cd8}.csv
#   03_results/03_pseudobulk/tables/pseudobulk_counts.csv        donor x condition x label, Ensembl-keyed
#   03_results/03_pseudobulk/tables/pseudobulk_coldata.csv
#   03_results/03_pseudobulk/tables/gene_symbols.csv             the Ensembl->HGNC seam map
#   03_results/05_scoring/tables/gsea_pseudobulk_{tag}.csv       the published NES (gate only)
#   00_data/references/{msigdb_hallmark,temp_hsr_lens,tcr_activation_lens}/*.txt
#                                                                frozen curated lists
#   01_modules/RNAseq-toolkit/data/references/mitocarta3.0/processed/Homo_sapiens/mito_mitopathways.rds
#   ../mouse_anchor/03_results/human_projection/signatures/*/*_up.txt
#   ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv
#   ../sting_positive_control/03_results/06_reference_axis/signatures/{sting_specific,ifn_only}_up.txt
#
# Outputs — 03_results/14_unbiased_enrichment/tables/:
#   geneset_manifest.csv            every collection, its size before/after filtering
#   ranked_list_keycheck.csv        the symbol-vs-Ensembl guard, per population
#   geneset_overlap.csv             per population x database gene overlap counts
#   geneset_alias_map.csv           every set scored twice, which copy the pooled family
#                                   kept, and both copies' statistics side by side
#   wt_heat_up_reproduction.csv     the gate: published NES vs this run's NES
#   gsea_<population>_<database>.csv  per-database results incl. leading edge
#   gsea_all.csv                    tidy sweep with padj_pooled + pooled family size,
#                                   one row per (population, pathway_id). Two
#                                   denominators, deliberately named apart:
#                                   n_sets_scored_in_db is what the per-database `padj`
#                                   was corrected over (alias copies included),
#                                   n_tests_in_db is what that database contributes to
#                                   the pooled family, n_tests_pooled the whole family.
#   gsea_pooled_summary_by_db.csv   per database, significant before/after pooling
#   runsum_interactive_<population>_<set>.csv  running-sum substrate (stage-05 schema),
#                                   for every mouse-derived arm, every set named in
#                                   unbiased_enrichment.runsum_always, and the top-N
#                                   curated per population; runsum_interactive_index.csv
#                                   flags which of the two reasons emitted each curve
#   progeny_activity.csv            PROGENy MLM on the moderated-t contrast statistics
#   progeny_donor_activity.csv      PROGENy MLM per donor-pseudobulk sample
#   progeny_sf_vs_pb.csv            donor-paired SF-vs-PB test per population x pathway
#
# Objects (checkpoints, not deliverables) — 03_results/objects/:
#   14_genesets.rds                 all collections as fgsea `pathways` lists
#   14_gsea/<tag>__<database>.rds   one clusterProfiler gseaResult per cell of the sweep
#   14_progeny.rds                  raw decoupleR MLM results
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/14_unbiased_enrichment.R
#
# COMPUTE ONLY — no ggplot/ggsave. Figures live in 14_unbiased_enrichment_viz.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(clusterProfiler)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, round_numeric_cols
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_processing/pathway_utils.R")  # list_to_term2gene

STAGE  <- "14_unbiased_enrichment"
SCRIPT <- "02_analysis/scripts/14_unbiased_enrichment.R"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# 0. CONFIG — every parameter read from analysis_config.yaml, none chosen here
# ============================================================================

CFG  <- FIG_CFG
UE   <- CFG$unbiased_enrichment
if (is.null(UE))
  stop("[14] analysis_config.yaml has no `unbiased_enrichment:` block — add it before running.")

THR      <- CFG$thresholds
MINSZ    <- as.integer(THR$gsea_min_size %||% 5L)
MAXSZ    <- as.integer(THR$gsea_max_size %||% 500L)
SEED     <- as.integer(THR$gsea_seed     %||% 123L)
NPERM    <- as.integer(THR$gsea_nperm    %||% 100000L)
FDR      <- as.numeric(THR$gsea_fdr      %||% 0.05)
POOL_M   <- UE$padj_pooled_method %||% "BH"
SPECIES  <- CFG$project$species %||% "Homo sapiens"
RUNSUM_N <- as.integer(UE$runsum_top_curated %||% 5L)
RUNSUM_ALWAYS <- as.character(unlist(UE$runsum_always %||% character(0)))
NES_TOL  <- as.numeric(UE$reproduction_check$nes_tolerance %||% 0.01)

# Sorted populations: display name -> ranked-list tag -> contrast label. The contrast
# label matches the published stage-05 tables exactly so the running-sum substrates
# this stage emits are interchangeable with the ones already on disk.
POPS <- c(Treg = "treg", Tcon = "tcon", CD8 = "cd8")
contrast_label <- function(pop) sprintf("SF_vs_PB_%s", pop)

# The file-backed collections: frozen one-symbol-per-line lists this repo owns, as
# opposed to the collections fetched from msigdbr or read from the CollecTRI network.
# Named once because three separate rules key off the distinction — they are exempt
# from the nominal size filter, they are the only collections eligible to be demoted
# to a pooling alias (section 2b), and the manifest reports the exemption.
#
# ORDER IS LOAD-BEARING. It sets the declaration order these collections enter
# COLLECTIONS in, and section 2b keys canonicality off declaration order. `hsr_lens`
# comes after `project_frozen` because project_frozen already carries HSR_core as its
# seventh file: putting hsr_lens later leaves HSR_core canonical where it already was,
# so adding the collection does not move that set's pooled adjusted p.
FILE_BACKED_KEYS <- c("mouse_projection", "project_frozen", "sting_axes",
                      "hsr_lens", "tcr_activation")
FILE_BACKED_DBS  <- unlist(lapply(FILE_BACKED_KEYS,
                                  function(k) UE[[k]]$name %||% character(0)))

RESULTS   <- CFG$paths$results %||% "03_results/"
DIR_OBJ   <- CFG$paths$objects %||% "03_results/objects/"
TBL       <- file.path(RESULTS, STAGE, CFG$paths$stage_tables_subdir %||% "tables")
DIR_GSEA  <- file.path(DIR_OBJ, "14_gsea")
for (d in c(TBL, file.path(TBL, "_overview"), DIR_OBJ, DIR_GSEA))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

RANKED_DIR <- file.path(RESULTS, "03_pseudobulk",
                        CFG$paths$stage_tables_subdir %||% "tables")

message("=================================================================")
message("14_unbiased_enrichment — unbiased sweep of the JIA niche contrast")
message(sprintf("  species=%s  size=[%d,%d]  seed=%d  nperm=%d  pooled=%s",
                SPECIES, MINSZ, MAXSZ, SEED, NPERM, POOL_M))
message("=================================================================")

# ============================================================================
# 1. RANKED LISTS + the key check that stops a silent biological-null
# ============================================================================

#' Read one 2-column signed-statistic ranked list into a named, sorted vector.
#'
#' Deduplication keeps the FIRST occurrence of a symbol, which is the most extreme
#' |t| because the upstream R seam already collapsed to one row per symbol on max
#' |t| and wrote the list sorted. Identical to helpers/fgsea_prerank.R so the two
#' stages cannot drift.
read_ranked <- function(path) {
  if (!file.exists(path)) stop("[14] ranked list not found: ", path)
  d <- utils::read.table(path, sep = "\t", header = FALSE,
                         col.names = c("symbol", "stat"), stringsAsFactors = FALSE)
  d <- d[!is.na(d$stat) & nzchar(d$symbol), , drop = FALSE]
  d <- d[!duplicated(d$symbol), , drop = FALSE]
  v <- stats::setNames(d$stat, d$symbol)
  sort(v, decreasing = TRUE)
}

#' Confirm a ranked list is keyed by HGNC symbol, not Ensembl id.
#'
#' The failure this exists for is not a crash. An Ensembl-keyed list intersects
#' every reference collection at ~zero and fgsea returns empty/NA rows, which reads
#' downstream as "this contrast contains nothing" rather than as a broken join. So
#' the check is a hard stop with the diagnosis attached, and its numbers are
#' published in ranked_list_keycheck.csv rather than only logged.
keycheck_row <- function(pop, tag, ranked) {
  g <- names(ranked)
  frac_ens <- mean(grepl("^ENS[A-Z]*G[0-9]{6,}$", g))
  if (frac_ens > 0.5)
    stop(sprintf(paste0("[14] ranked_%s.tsv is keyed by Ensembl id (%.0f%% of keys look like ",
                        "ENSG..., e.g. %s). Every reference gene set matches on HGNC symbol, so ",
                        "fgsea would return empty/NA rows and the sweep would look like a ",
                        "biological null. Re-run 03b_pseudobulk_de.R with the gene_symbols.csv ",
                        "join before this stage."),
                 tag, 100 * frac_ens, g[1]))
  tibble::tibble(population = pop, tag = tag, n_ranked = length(ranked),
                 frac_keys_ensembl_like = frac_ens,
                 first_key = g[1], last_key = g[length(g)],
                 stat_min = min(ranked), stat_max = max(ranked))
}

RANKED <- lapply(POPS, function(tag) read_ranked(file.path(RANKED_DIR, sprintf("ranked_%s.tsv", tag))))
names(RANKED) <- names(POPS)

keycheck <- dplyr::bind_rows(lapply(names(POPS), function(pop)
  keycheck_row(pop, POPS[[pop]], RANKED[[pop]])))
readr::write_csv(round_numeric_cols(keycheck), file.path(TBL, "ranked_list_keycheck.csv"))
for (i in seq_len(nrow(keycheck)))
  message(sprintf("[1] %-4s ranked list: %d symbols, top=%s, t in [%.2f, %.2f]",
                  keycheck$population[i], keycheck$n_ranked[i], keycheck$first_key[i],
                  keycheck$stat_min[i], keycheck$stat_max[i]))

# ============================================================================
# 2. GENE-SET COLLECTIONS — build once, cache, publish a manifest
# ============================================================================

#' Size-filter an fgsea `pathways` list on NOMINAL set size.
#'
#' clusterProfiler re-applies the same bounds to the EFFECTIVE size (after
#' intersecting the ranked list), so this pass only trims sets that could never
#' qualify. Both numbers are reported, because a set dropped for effective size is
#' a coverage fact about the data and a set dropped for nominal size is a curation
#' fact about the collection.
filter_by_size <- function(gsets, min_sz = MINSZ, max_sz = MAXSZ) {
  sz <- vapply(gsets, length, integer(1))
  gsets[sz >= min_sz & sz <= max_sz]
}

#' Fetch one MSigDB collection as a named list of HGNC-symbol vectors.
#'
#' msigdbr 26.1.0 renamed the arguments category/subcategory to
#' collection/subcollection and renamed the subcollection CP:KEGG to
#' CP:KEGG_LEGACY. Both are handled by detection and retry so the CONFIG can keep
#' the mouse-anchor spelling and the two compartments stay comparable at a glance.
load_msigdb <- function(category, subcategory = "", species = SPECIES) {
  has_new <- "collection" %in% names(formals(msigdbr::msigdbr))
  fetch <- function(sub_use) {
    if (has_new)
      msigdbr::msigdbr(species = species, collection = category,
                       subcollection = if (nzchar(sub_use)) sub_use else NULL)
    else
      msigdbr::msigdbr(species = species, category = category,
                       subcategory = if (nzchar(sub_use)) sub_use else NULL)
  }
  used <- subcategory
  df <- tryCatch(fetch(subcategory), error = function(e) {
    if (nzchar(subcategory) && grepl("[Uu]nknown subcollection", conditionMessage(e))) {
      alt <- sub("CP:KEGG$", "CP:KEGG_LEGACY", subcategory)
      if (alt != subcategory) {
        message(sprintf("    [msigdb] %s/%s -> retrying subcollection '%s' (msigdbr 26 rename)",
                        category, subcategory, alt))
        used <<- alt
        return(fetch(alt))
      }
    }
    stop("[14] MSigDB fetch failed for ", category, "/", subcategory, ": ",
         conditionMessage(e))
  })
  if (is.null(df) || nrow(df) == 0)
    stop("[14] MSigDB ", category, "/", subcategory, " returned 0 rows.")
  sym <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else "human_gene_symbol"
  out <- lapply(split(as.character(df[[sym]]), as.character(df$gs_name)), unique)
  # The RESOLVED subcollection and the release travel with the sets, because
  # CP:KEGG_LEGACY is a materially different collection from the CP:KEGG the config
  # names, and a provenance record that hid the substitution would be wrong.
  attr(out, "subcollection_used") <- used
  attr(out, "db_version") <- if ("db_version" %in% colnames(df))
    paste(unique(as.character(df$db_version)), collapse = ",") else NA_character_
  out
}

#' Read a one-symbol-per-line frozen gene list; the file stem is the set name.
read_set_file <- function(path) {
  if (!file.exists(path)) stop("[14] frozen gene list not found: ", path)
  g <- trimws(readLines(path, warn = FALSE))
  g <- unique(g[nzchar(g)])
  if (length(g) == 0) stop("[14] frozen gene list is empty: ", path)
  g
}

#' CollecTRI regulons -> one UNSIGNED gene set per TF.
#'
#' Unsigned on purpose: pooling activating and repressing targets makes the set the
#' TF's transcriptional neighbourhood, which is what an MSigDB TFT set is and what a
#' ranked-list enrichment statistic can speak to. The SIGNED use of this same
#' network is decoupleR ULM, a different method with a different readout; running it
#' here would silently mix the two.
load_collectri <- function(spec) {
  if (!file.exists(spec$path)) stop("[14] CollecTRI network not found: ", spec$path)
  d <- readr::read_csv(spec$path, show_col_types = FALSE, progress = FALSE)
  sc <- spec$source_col %||% "source"; tc <- spec$target_col %||% "target"
  if (!all(c(sc, tc) %in% colnames(d)))
    stop("[14] CollecTRI table lacks columns ", sc, "/", tc)
  lapply(split(as.character(d[[tc]]), as.character(d[[sc]])), unique)
}

#' Read one toolkit reference-database RDS as a named list of gene-symbol vectors.
#'
#' The toolkit ships these pre-converted per species as a list carrying a T2G frame
#' (gs_name, gene_symbol), so the species is a property of the FILE PATH and not of
#' anything this function can see. That is the whole failure mode worth guarding: a
#' Mus_musculus path silently substituted here would deliver Title-cased mouse symbols,
#' intersect this compartment's HGNC-keyed ranked lists at approximately zero, and fgsea
#' would report empty rather than wrong. So the loader asserts the symbols look human
#' and stops with the path in the message when they do not.
load_reference_rds <- function(spec) {
  if (!file.exists(spec$path))
    stop("[14] reference database RDS not found: ", spec$path)
  obj <- readRDS(spec$path)
  if (!is.list(obj) || !"T2G" %in% names(obj))
    stop("[14] ", spec$path, " is not a toolkit reference database (no T2G element).")
  t2g <- obj$T2G
  if (ncol(t2g) < 2)
    stop("[14] ", spec$path, " T2G needs a term column and a gene column.")
  sets <- lapply(split(as.character(t2g[[2]]), as.character(t2g[[1]])), unique)
  sets <- lapply(sets, function(g) g[!is.na(g) & nzchar(g)])
  sets <- sets[lengths(sets) > 0]
  if (length(sets) == 0) stop("[14] ", spec$path, " yielded no non-empty gene set.")
  all_sym  <- unique(unlist(sets, use.names = FALSE))
  frac_upr <- mean(all_sym == toupper(all_sym))
  if (frac_upr < 0.9)
    stop(sprintf(paste0("[14] %s carries symbols that are mostly not uppercase (%.0f%% are, ",
                        "e.g. %s), so this is the mouse-symbol build. Every reference set in ",
                        "this compartment must match the HGNC-keyed ranked lists, and a mouse ",
                        "build would intersect them at approximately zero — which fgsea reports ",
                        "as an empty result rather than as an error. Point `path` at the ",
                        "processed/Homo_sapiens build."),
                 spec$path, 100 * frac_upr,
                 paste(utils::head(all_sym[all_sym != toupper(all_sym)], 4), collapse = ", ")))
  attr(sets, "db_source") <- sprintf("%s (%s)",
                                     as.character(obj$source %||% basename(spec$path)),
                                     spec$path)
  sets
}

# Every collection the CONFIG asks for, in declaration order. Computed before the cache
# is consulted, because the cache is only valid if it holds exactly this set of names.
CONFIGURED_DBS <- c(
  vapply(UE$msigdb, function(m) as.character(m$name), character(1)),
  UE$tf_network$name %||% character(0),
  unlist(lapply(UE$custom_rds %||% list(), function(s) as.character(s$name))),
  unlist(lapply(FILE_BACKED_KEYS, function(k) UE[[k]]$name %||% character(0))))

genesets_path <- file.path(DIR_OBJ, "14_genesets.rds")
## THE CACHE IS VALIDATED AGAINST THE CONFIG, not merely tested for existence. Adding a
## database to analysis_config.yaml with a cache already on disk used to load the old
## collections and sweep those instead, and every table would then report the previous
## family under the new config with nothing to say so. The mismatch is a rebuild, not an
## error, because rebuilding is cheap and correct.
or_none <- function(x) if (length(x) == 0) "none" else paste(x, collapse = ", ")
cache_dbs <- if (file.exists(genesets_path)) names(readRDS(genesets_path)) else character(0)
cache_ok  <- file.exists(genesets_path) && setequal(cache_dbs, CONFIGURED_DBS)
if (file.exists(genesets_path) && !cache_ok)
  message(sprintf(paste0("[2] gene-set cache is stale: it holds %d collection(s) and the config ",
                         "asks for %d (added: %s; removed: %s). Rebuilding."),
                  length(cache_dbs), length(CONFIGURED_DBS),
                  or_none(setdiff(CONFIGURED_DBS, cache_dbs)),
                  or_none(setdiff(cache_dbs, CONFIGURED_DBS))))

if (cache_ok) {
  message("[2] gene-set collections: cache hit -> ", genesets_path)
  COLLECTIONS <- readRDS(genesets_path)
} else {
  message("[2] building gene-set collections ...")
  COLLECTIONS <- list()

  for (m in UE$msigdb) {
    raw  <- load_msigdb(m$category, m$subcategory %||% "")
    used <- attr(raw, "subcollection_used") %||% ""
    src  <- sprintf("msigdbr %s / MSigDB %s, collection %s%s",
                    as.character(utils::packageVersion("msigdbr")),
                    attr(raw, "db_version") %||% "unknown", m$category,
                    if (nzchar(used)) paste0("/", used) else "")
    if (!identical(used, m$subcategory %||% ""))
      src <- paste0(src, sprintf(" (config asked for %s)", m$subcategory))
    COLLECTIONS[[m$name]] <- list(sets = filter_by_size(raw), n_raw = length(raw),
                                  source = src)
    message(sprintf("  %-14s %5d raw -> %5d kept   [%s]", m$name, length(raw),
                    length(COLLECTIONS[[m$name]]$sets), src))
  }

  raw <- load_collectri(UE$tf_network)
  COLLECTIONS[[UE$tf_network$name]] <- list(sets = filter_by_size(raw), n_raw = length(raw),
                                            source = sprintf("CollecTRI regulons, unsigned (%s)",
                                                             UE$tf_network$path))
  message(sprintf("  %-14s %5d raw -> %5d kept", UE$tf_network$name, length(raw),
                  length(COLLECTIONS[[UE$tf_network$name]]$sets)))

  # Toolkit reference databases shipped pre-converted per species. Nominal-size filtered
  # like the msigdbr collections, because they are curated catalogues of comparable width
  # rather than the handful of sets a reader came for.
  for (spec in UE$custom_rds %||% list()) {
    raw <- load_reference_rds(spec)
    COLLECTIONS[[spec$name]] <- list(sets = filter_by_size(raw), n_raw = length(raw),
                                     source = attr(raw, "db_source"))
    message(sprintf("  %-14s %5d raw -> %5d kept   [%s]", spec$name, length(raw),
                    length(COLLECTIONS[[spec$name]]$sets), attr(raw, "db_source")))
  }

  # File-backed collections: mouse-derived UP arms, this compartment's frozen curated
  # lists, and the frozen SAVI axes. These are NOT nominal-size filtered — they are
  # the sets the reader came for, and a set below the floor must be reported as
  # untestable WITH its size rather than vanish. clusterProfiler still declines to
  # score one, and its absence from the results is recorded in geneset_overlap.csv.
  for (key in FILE_BACKED_KEYS) {
    spec <- UE[[key]]
    if (is.null(spec)) next
    sets <- stats::setNames(
      lapply(spec$files, read_set_file),
      vapply(spec$files, function(p) sub("\\.txt$", "", basename(p)), character(1)))
    COLLECTIONS[[spec$name]] <- list(
      sets = sets, n_raw = length(sets),
      source = sprintf("frozen gene lists: %s", paste(spec$files, collapse = "; ")))
    message(sprintf("  %-14s %5d sets (%s)", spec$name, length(sets),
                    paste(sprintf("%s=%d", names(sets), lengths(sets)), collapse = " ")))
  }

  saveRDS(COLLECTIONS, genesets_path)
  message("  cached -> ", genesets_path)
}

DBS <- names(COLLECTIONS)
GATE_DB <- UE$mouse_projection$name %||% "mouse_projection"
if (!GATE_DB %in% DBS) stop("[14] gate database '", GATE_DB, "' missing from the collections.")
# Cheap collections first so a failure is cheap; GO_BP (the largest) last.
DB_ORDER <- c(GATE_DB, setdiff(DBS, GATE_DB))
DB_ORDER <- c(DB_ORDER[1], DB_ORDER[-1][order(vapply(COLLECTIONS[DB_ORDER[-1]],
                                                     function(x) length(x$sets), integer(1)))])

# ============================================================================
# 2b. ALIAS RESOLUTION — one set, one hypothesis in the pooled family
# ============================================================================
## `project_frozen` re-pins six MSigDB Hallmark sets to files, and the `Hallmark`
## collection fetches those same six from msigdbr. The gene content is identical, so
## before this resolution existed the pooled Benjamini-Hochberg family contained six
## exact duplicate hypotheses per population and every rank denominator was six too
## high. The frozen files are NOT removed from the config: they are what
## freeze_hallmark_sets.R produced and what the decomposition and purge stages actually
## consume, so a manifest that stopped listing them would hide a real dependency. The
## collision is resolved instead — each copy is still scored and still appears in its
## own per-database table, but only the canonical copy enters the pooled family.
##
## Detection is by ID COLLISION PLUS GENE-CONTENT IDENTITY, never by a hardcoded list of
## set names, so a seventh frozen Hallmark set added to the config tomorrow is
## deduplicated by construction without editing this file. Two situations are hard stops
## rather than silent resolutions:
##   - same id, DIFFERENT genes. The id no longer names one hypothesis, and choosing a
##     winner would silently discard a real set. A human has to look.
##   - a collection that is NOT file-backed would be the one demoted. MSigDB is the
##     canonical home of a HALLMARK_/GOBP_/REACTOME_ id and a frozen file in this repo
##     is the re-pin of it, so only FILE_BACKED_DBS are alias-eligible.
##
## Precedence is the config DECLARATION order, not DB_ORDER. DB_ORDER is a run-cost
## heuristic that puts the smallest collection first; keying canonicality off it would
## make which copy survives depend on how many sets a collection happens to contain.

DECL_ORDER <- names(COLLECTIONS)
set_owners <- list()
for (db in DECL_ORDER)
  for (id in names(COLLECTIONS[[db]]$sets))
    set_owners[[id]] <- c(set_owners[[id]], db)

alias_rows <- list()
for (id in names(set_owners)[vapply(set_owners, length, integer(1)) > 1L]) {
  dbs   <- set_owners[[id]]
  genes <- lapply(dbs, function(d) sort(unique(COLLECTIONS[[d]]$sets[[id]])))
  if (!all(vapply(genes[-1], function(g) identical(g, genes[[1]]), logical(1))))
    stop(sprintf(paste0("[14] set id '%s' is present in %s with DIFFERENT gene content ",
                        "(%s). The id therefore does not name one hypothesis, and ",
                        "deduplicating it would silently discard a real set. Rename one of ",
                        "them in analysis_config.yaml, or decide which is intended."),
                 id, paste(dbs, collapse = " and "),
                 paste(sprintf("%s n=%d", dbs, lengths(genes)), collapse = ", ")))
  canon   <- dbs[1]
  aliases <- dbs[-1]
  demoted <- setdiff(aliases, FILE_BACKED_DBS)
  if (length(demoted))
    stop(sprintf(paste0("[14] resolving set id '%s' would demote the curated collection(s) ",
                        "%s to an alias of %s. Only the file-backed collections (%s) are ",
                        "alias-eligible, because MSigDB is the canonical home of a curated ",
                        "set id and a frozen file here is the re-pin of it. Check the order ",
                        "collections are declared in analysis_config.yaml."),
                 id, paste(demoted, collapse = ", "), canon,
                 paste(FILE_BACKED_DBS, collapse = ", ")))
  for (a in aliases)
    alias_rows[[length(alias_rows) + 1L]] <- tibble::tibble(
      pathway_id = id, alias_database = a, canonical_database = canon,
      n_genes = length(genes[[1]]), gene_content_identical = TRUE,
      resolution = "scored in both databases; enters the pooled family under canonical_database only")
}
ALIAS <- dplyr::bind_rows(alias_rows)

#' Is this (database, set) copy an alias that the pooled family must exclude?
#'
#' Keyed on the resolved pair, so the pooling step excludes by construction from
#' section 2b's structural detection rather than by re-matching set names downstream.
ALIAS_KEY <- if (nrow(ALIAS) > 0) paste(ALIAS$alias_database, ALIAS$pathway_id) else character(0)
is_pooling_alias <- function(db, id) paste(db, id) %in% ALIAS_KEY

if (nrow(ALIAS) == 0) {
  message("[2b] alias resolution: no set id appears in more than one collection.")
} else {
  message(sprintf("[2b] alias resolution: %d set(s) appear in more than one collection; each is scored in both and pooled once.",
                  nrow(ALIAS)))
  for (i in seq_len(nrow(ALIAS)))
    message(sprintf("     %-42s pooled under %-14s alias copy in %-14s (%d genes, identical)",
                    ALIAS$pathway_id[i], ALIAS$canonical_database[i],
                    ALIAS$alias_database[i], ALIAS$n_genes[i]))
}

manifest <- dplyr::bind_rows(lapply(DB_ORDER, function(db) {
  x <- COLLECTIONS[[db]]
  n_alias <- if (nrow(ALIAS) > 0) sum(ALIAS$alias_database == db) else 0L
  tibble::tibble(database = db, source = x$source, n_sets_in_source = x$n_raw,
                 n_sets_after_nominal_size_filter = length(x$sets),
                 nominal_size_filter_applied = !db %in% FILE_BACKED_DBS,
                 n_sets_aliased_out_of_pooling = n_alias,
                 n_sets_offered_for_pooling = length(x$sets) - n_alias,
                 gsea_min_size = MINSZ, gsea_max_size = MAXSZ,
                 min_set_size = min(lengths(x$sets)), max_set_size = max(lengths(x$sets)))
}))
readr::write_csv(round_numeric_cols(manifest), file.path(TBL, "geneset_manifest.csv"))
message(sprintf("[2] manifest -> %s (%d databases, %d sets total)",
                file.path(TBL, "geneset_manifest.csv"), nrow(manifest),
                sum(manifest$n_sets_after_nominal_size_filter)))

# ============================================================================
# 3. OVERLAP — the published evidence that the symbol join actually joined
# ============================================================================

overlap <- dplyr::bind_rows(lapply(names(POPS), function(pop) {
  g <- names(RANKED[[pop]])
  dplyr::bind_rows(lapply(DB_ORDER, function(db) {
    u  <- unique(unlist(COLLECTIONS[[db]]$sets, use.names = FALSE))
    ov <- length(intersect(u, g))
    tibble::tibble(population = pop, database = db, n_ranked = length(g),
                   n_unique_set_genes = length(u), n_overlap = ov,
                   frac_of_set_genes = ov / length(u), frac_of_ranked = ov / length(g))
  }))
}))
readr::write_csv(round_numeric_cols(overlap), file.path(TBL, "geneset_overlap.csv"))

# Hard guard on the largest curated collection actually present. A handful of
# matching symbols means the join failed; fgsea would not say so.
gbp <- overlap[overlap$database == "GO_BP", ]
if (nrow(gbp) > 0 && min(gbp$n_overlap) < 2000)
  stop(sprintf(paste0("[14] GO_BP overlaps the ranked lists at only %d genes (minimum across ",
                      "populations). Expected many thousands. This is the Ensembl-versus-HGNC ",
                      "silent failure, not a biological null — check the ranked-list keys."),
               min(gbp$n_overlap)))
for (pop in names(POPS)) {
  o <- overlap[overlap$population == pop, ]
  message(sprintf("[3] %-4s overlap: %s", pop,
                  paste(sprintf("%s=%d", o$database, o$n_overlap), collapse = " ")))
}

# ============================================================================
# 4. fgsea ENGINE — identical settings to helpers/fgsea_prerank.R
# ============================================================================

#' Run clusterProfiler::GSEA(by="fgsea") for one (population, database) cell.
#'
#' Every engine setting is the one helpers/fgsea_prerank.R used for the published
#' stage-05 run — exponent 1, eps 0 (exact multilevel p-values), pvalueCutoff 1 (keep
#' every set; FDR filtering is a reporting choice, not a run-time one), BH within the
#' call, seeded. Deviating on any of these is what would break the reproduction gate.
run_cell <- function(ranked, sets, db, pop) {
  cache <- file.path(DIR_GSEA, sprintf("%s__%s.rds", POPS[[pop]], db))
  if (file.exists(cache)) {
    message(sprintf("  [gsea] %-5s x %-14s cache hit", pop, db))
    return(readRDS(cache))
  }
  message(sprintf("  [gsea] %-5s x %-14s %d sets ...", pop, db, length(sets)))
  t0 <- Sys.time()
  set.seed(SEED)
  g <- clusterProfiler::GSEA(
    geneList      = ranked,
    TERM2GENE     = list_to_term2gene(sets),
    by            = "fgsea",
    exponent      = 1,
    eps           = 0,
    minGSSize     = MINSZ,
    maxGSSize     = MAXSZ,
    nPermSimple   = NPERM,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    seed          = TRUE,
    verbose       = FALSE
  )
  message(sprintf("         -> %d scored in %.0f s", nrow(g@result),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(g, cache)
  g
}

#' Project one gseaResult onto the compartment's master_gsea_table schema plus the
#' sweep columns. `core_enrichment` (the leading edge) is kept here, in the
#' per-database file; the tidy cross-database table carries only its SIZE so it stays
#' readable, and points back to these files for the gene names.
#'
#' `is_pooled_alias` marks a row whose set is scored here but pooled under another
#' collection (section 2b). It is carried in the per-database files, which is where a
#' reader meets those rows; gsea_all.csv drops them, so every id appears once there.
tidy_cell <- function(g, db, pop) {
  r <- g@result
  if (nrow(r) == 0) return(NULL)
  tibble::tibble(
    pathway_id         = r$ID,
    pathway_name       = r$Description,
    database           = db,
    population         = pop,
    contrast           = contrast_label(pop),
    nes                = r$NES,
    enrichment_score   = r$enrichmentScore,
    pvalue             = r$pvalue,
    padj               = r$p.adjust,
    set_size           = r$setSize,
    leading_edge_size  = vapply(strsplit(r$core_enrichment, "/", fixed = TRUE),
                                function(x) length(x[nzchar(x)]), integer(1)),
    core_enrichment    = r$core_enrichment,
    direction          = ifelse(r$NES > 0, "up", "down"),
    is_pooled_alias    = vapply(r$ID, function(i) is_pooling_alias(db, i), logical(1),
                                USE.NAMES = FALSE)
  ) |> dplyr::arrange(pvalue, dplyr::desc(abs(nes)))
}

# ============================================================================
# 5. THE CORRECTNESS GATE — run before the expensive collections
# ============================================================================
## WT_heat_up is scored as an ordinary member of the mouse_projection collection and
## its NES compared against the published stage-05 value off the SAME ranked list.
## A mismatch means this stage's set prep, ranking or engine settings differ from the
## published pipeline; the sweep would then not be comparable to the targeted result
## it exists to calibrate, so the run stops here. It is never retuned to agree.

message("[5] reproduction gate: WT_heat_up against the published stage-05 NES ...")
GATE_SIG   <- UE$reproduction_check$signature %||% "WT_heat_up"
gate_rows  <- list()
gsea_cells <- list()   # (pop, db) -> gseaResult

for (pop in names(POPS)) {
  g <- run_cell(RANKED[[pop]], COLLECTIONS[[GATE_DB]]$sets, GATE_DB, pop)
  gsea_cells[[paste(pop, GATE_DB)]] <- g

  pub_path <- file.path(RESULTS, "05_scoring", CFG$paths$stage_tables_subdir %||% "tables",
                        sprintf("gsea_pseudobulk_%s.csv", POPS[[pop]]))
  if (!file.exists(pub_path))
    stop("[14] published stage-05 table absent, cannot verify the gate: ", pub_path)
  pub <- readr::read_csv(pub_path, show_col_types = FALSE, progress = FALSE)
  pr  <- pub[pub$pathway_id == GATE_SIG, ]
  hr  <- g@result[g@result$ID == GATE_SIG, ]
  if (nrow(pr) != 1 || nrow(hr) != 1)
    stop(sprintf("[14] gate: %s found %d time(s) published and %d time(s) here for %s.",
                 GATE_SIG, nrow(pr), nrow(hr), pop))

  gate_rows[[pop]] <- tibble::tibble(
    population = pop, signature = GATE_SIG,
    nes_published = pr$nes[1], nes_this_stage = hr$NES[1],
    abs_nes_difference = abs(pr$nes[1] - hr$NES[1]),
    set_size_published = pr$set_size[1], set_size_this_stage = hr$setSize[1],
    pvalue_published = pr$pvalue[1], pvalue_this_stage = hr$pvalue[1],
    nes_tolerance = NES_TOL,
    reproduced = abs(pr$nes[1] - hr$NES[1]) <= NES_TOL &&
                 pr$set_size[1] == hr$setSize[1])
  message(sprintf("    %-4s published NES %+.6f  this stage %+.6f  |diff| %.2e  setSize %d vs %d  -> %s",
                  pop, pr$nes[1], hr$NES[1], abs(pr$nes[1] - hr$NES[1]),
                  pr$set_size[1], hr$setSize[1],
                  ifelse(gate_rows[[pop]]$reproduced, "REPRODUCED", "MISMATCH")))
}

gate <- dplyr::bind_rows(gate_rows)
readr::write_csv(round_numeric_cols(gate), file.path(TBL, "wt_heat_up_reproduction.csv"))
if (!all(gate$reproduced))
  stop(sprintf(paste0("[14] REPRODUCTION GATE FAILED for %s. Published NES %s versus this ",
                      "stage %s (tolerance %g). The set prep, ranking or fgsea parameters ",
                      "differ from the published pipeline, so the whole sweep would be on a ",
                      "different footing. Diagnose the difference; do NOT tune parameters ",
                      "until it agrees."),
               paste(gate$population[!gate$reproduced], collapse = ", "),
               paste(sprintf("%.6f", gate$nes_published[!gate$reproduced]), collapse = "/"),
               paste(sprintf("%.6f", gate$nes_this_stage[!gate$reproduced]), collapse = "/"),
               NES_TOL))
message("[5] gate passed in every population.")

# ============================================================================
# 6. THE SWEEP — remaining databases, per population
# ============================================================================

message("[6] sweeping the remaining databases ...")
for (pop in names(POPS)) {
  for (db in setdiff(DB_ORDER, GATE_DB))
    gsea_cells[[paste(pop, db)]] <- run_cell(RANKED[[pop]], COLLECTIONS[[db]]$sets, db, pop)
}

# ============================================================================
# 7. PER-DATABASE TABLES + pooled BH across the whole family per population
# ============================================================================
## Two multiplicity corrections are reported side by side, because they answer
## different questions. `padj` is BH WITHIN one database, which is what a
## single-collection run would report and what makes a row comparable to a published
## per-collection result. `padj_pooled` is BH across EVERY test asked of that
## population's ranked list, which is the honest correction for a sweep this wide:
## fourteen databases interrogated with one ranking is one family of hypotheses, not
## fourteen independent studies. A row that survives only the per-database correction
## is a row whose significance depends on not counting the rest of the sweep.
##
## The pooled family is DEDUPLICATED FIRST. Section 2b resolved which copy of a set that
## appears in two collections is canonical; the alias copies are dropped here before
## p.adjust, so the family contains each hypothesis once and `n_tests_pooled` is the
## rank denominator a reader can quote. Per-database `padj` is computed on the full
## per-database table (that is what a standalone single-collection run reports, aliases
## and all), so the two corrections keep answering their two different questions.

message("[7] assembling tables and applying pooled multiplicity correction ...")
tidy_all <- list()
for (pop in names(POPS)) {
  for (db in DB_ORDER) {
    tt <- tidy_cell(gsea_cells[[paste(pop, db)]], db, pop)
    if (is.null(tt)) {
      warning(sprintf("[14] %s x %s scored 0 sets — written as an empty table.", pop, db),
              call. = FALSE)
      next
    }
    readr::write_csv(round_numeric_cols(tt),
                     file.path(TBL, sprintf("gsea_%s_%s.csv", POPS[[pop]], db)))
    tidy_all[[paste(pop, db)]] <- tt
  }
}

scored_all <- dplyr::bind_rows(tidy_all)
if (nrow(scored_all) == 0) stop("[14] the sweep produced no rows at all — refusing to continue.")

# TWO DENOMINATORS, NAMED APART, because mixing them prints a wrong number that reads
# as a real result. `n_sets_scored_in_db` counts every set scored in that database
# including alias copies — it is what the per-database `padj` in gsea_<pop>_<db>.csv was
# corrected over. `n_tests_in_db` counts only what that database contributes to the
# pooled family, so it is the denominator that belongs beside any pooled count. The
# overview figure annotates "<pooled-significant> of <n_tests_in_db>"; pairing a
# deduplicated numerator with the pre-dedup denominator made project_frozen read "0 of
# 7" when only one of its seven sets is in the pooled family at all.
scored_all <- scored_all |>
  dplyr::group_by(population, database) |>
  dplyr::mutate(n_sets_scored_in_db = dplyr::n()) |>
  dplyr::ungroup()

n_dropped <- sum(scored_all$is_pooled_alias)
sweep_df <- scored_all |>
  dplyr::filter(!is_pooled_alias) |>
  dplyr::group_by(population) |>
  dplyr::mutate(padj_pooled    = stats::p.adjust(pvalue, method = POOL_M),
                n_tests_pooled = dplyr::n()) |>
  dplyr::group_by(population, database) |>
  dplyr::mutate(n_tests_in_db  = dplyr::n()) |>
  dplyr::ungroup()
message(sprintf("  pooled family: %d alias row(s) dropped before BH; %d rows pooled",
                n_dropped, nrow(sweep_df)))

sweep_out <- sweep_df |>
  dplyr::select(population, contrast, database, pathway_id, pathway_name, direction,
                nes, enrichment_score, pvalue, padj, padj_pooled, set_size,
                leading_edge_size, n_sets_scored_in_db, n_tests_in_db, n_tests_pooled) |>
  dplyr::arrange(population, padj_pooled, dplyr::desc(abs(nes)))
readr::write_csv(round_numeric_cols(sweep_out), file.path(TBL, "gsea_all.csv"))
message(sprintf("  gsea_all.csv: %d rows (%s)", nrow(sweep_out),
                paste(sprintf("%s n=%d", names(POPS),
                              vapply(names(POPS), function(p) sum(sweep_df$population == p), integer(1))),
                      collapse = ", ")))

# The alias record, with BOTH copies' statistics side by side. Written here rather than
# in section 2b because the point of the table is that a reader can see the dropped copy
# agreed with the kept one to fgsea's permutation noise — a difference larger than that
# would mean the two copies are not the same test after all.
if (nrow(ALIAS) > 0) {
  alias_stat <- function(rows, db_col) {
    dplyr::select(rows, population, contrast, pathway_id,
                  !!paste0("nes_", db_col) := nes, !!paste0("pvalue_", db_col) := pvalue,
                  !!paste0("padj_in_", db_col) := padj, set_size)
  }
  alias_out <- ALIAS |>
    dplyr::cross_join(dplyr::distinct(scored_all, population, contrast)) |>
    dplyr::left_join(dplyr::select(scored_all, population, pathway_id, database,
                                  nes_kept = nes, pvalue_kept = pvalue, padj_kept = padj,
                                  set_size),
                     by = c("population", "pathway_id",
                            "canonical_database" = "database")) |>
    dplyr::left_join(dplyr::select(scored_all, population, pathway_id, database,
                                  nes_dropped = nes, pvalue_dropped = pvalue,
                                  padj_dropped = padj),
                     by = c("population", "pathway_id", "alias_database" = "database")) |>
    dplyr::mutate(abs_nes_difference = abs(nes_kept - nes_dropped),
                  kept_copy = canonical_database, dropped_copy = alias_database) |>
    dplyr::select(population, contrast, pathway_id, n_genes, gene_content_identical,
                  kept_copy, dropped_copy, set_size, nes_kept, nes_dropped,
                  abs_nes_difference, pvalue_kept, pvalue_dropped, padj_kept,
                  padj_dropped, resolution) |>
    dplyr::arrange(population, pathway_id)
  readr::write_csv(round_numeric_cols(alias_out), file.path(TBL, "geneset_alias_map.csv"))
  message(sprintf("  geneset_alias_map.csv: %d row(s); largest |NES| difference between copies %.2e",
                  nrow(alias_out), max(alias_out$abs_nes_difference, na.rm = TRUE)))
} else {
  readr::write_csv(tibble::tibble(population = character(), contrast = character(),
                                  pathway_id = character(), n_genes = integer(),
                                  gene_content_identical = logical(), kept_copy = character(),
                                  dropped_copy = character(), set_size = integer(),
                                  nes_kept = numeric(), nes_dropped = numeric(),
                                  abs_nes_difference = numeric(), pvalue_kept = numeric(),
                                  pvalue_dropped = numeric(), padj_kept = numeric(),
                                  padj_dropped = numeric(), resolution = character()),
                  file.path(TBL, "geneset_alias_map.csv"))
  message("  geneset_alias_map.csv: empty (no set appears in two collections).")
}

## Counted on the FULL scored table with `padj_pooled` joined back on, so
## `sig_per_database` is the count a standalone single-collection run would report
## (aliases included, matching that database's own CSV) while `sig_pooled` is counted on
## the deduplicated family. An alias row carries padj_pooled = NA by construction and so
## contributes to neither pooled count; `n_sets_aliased_out_of_pooling` says how many.
## Each significance count is divided by ITS OWN denominator — `sig_per_database` by
## `n_sets_scored_in_db`, `sig_pooled` by `n_tests_in_db` — so no row of this table
## pairs a count with a denominator it was not computed over.
scored_pooled <- scored_all |>
  dplyr::left_join(dplyr::select(sweep_df, population, database, pathway_id,
                                 padj_pooled, n_tests_pooled),
                   by = c("population", "database", "pathway_id"))

# The pooled family size is a property of the POPULATION, joined on rather than reduced
# out of each database's rows: a database whose every set were an alias would have no
# non-NA value to reduce, and a silent -Inf there is exactly the kind of quiet wrong
# number this stage exists to avoid.
pooled_family_n <- dplyr::distinct(sweep_df, population, n_tests_pooled)

pooled_by_db <- scored_pooled |>
  dplyr::group_by(population, contrast, database) |>
  dplyr::summarise(n_sets_scored_in_db = dplyr::n(),
                   n_sets_aliased_out_of_pooling = sum(is_pooled_alias),
                   n_tests_in_db      = sum(!is_pooled_alias),
                   sig_per_database   = sum(padj < FDR, na.rm = TRUE),
                   sig_pooled         = sum(padj_pooled < FDR, na.rm = TRUE),
                   # Counted over NON-ALIAS rows on both sides, so it means what it says:
                   # sets this database lost by facing the wider family. Counting alias
                   # rows as "lost" would blame pooling for six sets that are
                   # pooled-significant under their canonical collection.
                   sig_lost_to_pooling = sum(padj < FDR & !is_pooled_alias, na.rm = TRUE) -
                                         sum(padj_pooled < FDR, na.rm = TRUE),
                   sig_pooled_up      = sum(padj_pooled < FDR & nes > 0, na.rm = TRUE),
                   sig_pooled_down    = sum(padj_pooled < FDR & nes < 0, na.rm = TRUE),
                   pct_sig_per_database = round(100 * sum(padj < FDR, na.rm = TRUE) /
                                                dplyr::n(), 1),
                   pct_sig_pooled       = round(100 * sum(padj_pooled < FDR, na.rm = TRUE) /
                                                max(sum(!is_pooled_alias), 1L), 1),
                   top_pooled_up      = if (any(padj_pooled < FDR & nes > 0, na.rm = TRUE))
                     pathway_id[which.max(ifelse(padj_pooled < FDR & nes > 0, nes, -Inf))] else NA_character_,
                   top_pooled_up_nes  = if (any(padj_pooled < FDR & nes > 0, na.rm = TRUE))
                     max(nes[padj_pooled < FDR & nes > 0], na.rm = TRUE) else NA_real_,
                   min_pvalue         = min(pvalue, na.rm = TRUE),
                   n_at_min_pvalue    = sum(pvalue == min(pvalue, na.rm = TRUE), na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::left_join(pooled_family_n, by = "population") |>
  dplyr::relocate(n_tests_pooled, .after = n_tests_in_db) |>
  dplyr::arrange(population, dplyr::desc(sig_pooled))
readr::write_csv(round_numeric_cols(pooled_by_db),
                 file.path(TBL, "gsea_pooled_summary_by_db.csv"))

for (pop in names(POPS)) {
  s <- pooled_by_db[pooled_by_db$population == pop, ]
  message(sprintf("  %-4s pooled family = %d tests; significant %d per-database -> %d pooled at FDR<%.2g",
                  pop, s$n_tests_pooled[1], sum(s$sig_per_database), sum(s$sig_pooled), FDR))
}

# ============================================================================
# 8. RUNNING-SUM SUBSTRATE — schema-identical to the stage-05 interactive tables
# ============================================================================
## The downstream interactive running-sum comparison reads these, so the column set,
## order and names must match 05_scoring/tables/runsum_interactive_*.csv EXACTLY.
## The curve is recomputed with the DOSE weighted-KS formula off the fitted object's
## own geneList and exponent, so the emitted curve IS the plotted one; the
## leading_edge flag comes from the object's own core_enrichment.

#' Weighted running enrichment score for one set against one ranked vector.
runsum_table <- function(g, set_id, genes, pop) {
  gl <- g@geneList; n <- length(gl); gn <- names(gl)
  ex <- as.numeric(g@params[["exponent"]] %||% 1)
  hits <- gn %in% intersect(genes, gn)
  nh <- sum(hits)
  phit <- numeric(n); pmiss <- numeric(n)
  nr <- sum(abs(gl[hits])^ex)
  if (nr > 0) phit[hits] <- (abs(gl[hits])^ex) / nr
  pmiss[!hits] <- 1 / max(n - nh, 1L)
  core <- g@result[g@result$ID == set_id, "core_enrichment"]
  core <- if (length(core) == 1 && !is.na(core) && nzchar(core))
    strsplit(core, "/", fixed = TRUE)[[1]] else character(0)
  data.frame(rank = seq_len(n), gene = gn, stat = as.numeric(gl),
             running_es = cumsum(phit) - cumsum(pmiss), hit = hits,
             leading_edge = hits & (gn %in% core), gene_set = set_id,
             population = POPS[[pop]], contrast = contrast_label(pop),
             stringsAsFactors = FALSE)
}

message(sprintf("[8] emitting running-sum substrate (mouse up arms + %d named set(s) + top %d curated per population) ...",
                length(RUNSUM_ALWAYS), RUNSUM_N))
CURATED_DBS <- setdiff(DB_ORDER, GATE_DB)
RUNSUM_HALF <- max(as.integer(RUNSUM_N / 2), 1L)

## The canonical collection of a named set, resolved from the pooled table rather than by
## searching the collections: an id re-pinned into two collections has a row in both, and
## only the canonical one carries a pooled adjusted p (section 2b sets the alias copy's to
## NA). Looking the id up here therefore CANNOT emit two identical curves under two
## database names, which is the failure the alias resolution exists to prevent.
canonical_db_of <- function(id, pop) {
  hit <- sweep_df$database[sweep_df$population == pop & sweep_df$pathway_id == id &
                           !is.na(sweep_df$padj_pooled)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

runsum_index <- list()
for (pop in names(POPS)) {
  # Always: every mouse-derived up arm that was scorable.
  wanted <- lapply(names(COLLECTIONS[[GATE_DB]]$sets), function(s) list(db = GATE_DB, id = s))
  # Plus the config's named sets, in every population regardless of rank. These are the
  # comparators the compartment's question turns on, and a rank-based quota cannot
  # guarantee them: HALLMARK_HYPOXIA sits outside the top-N by pooled p in all three
  # populations, so the co-imposed niche stress had no curve to compare the mouse arm
  # against. Naming them here is a REPORTING choice and touches no statistic — the
  # trace is a deterministic walk of a ranking that is already fixed.
  for (id in RUNSUM_ALWAYS) {
    db <- canonical_db_of(id, pop)
    if (is.na(db)) {
      message(sprintf("    %s x %s: named in runsum_always but absent from the pooled family — no substrate",
                      pop, id))
      next
    }
    wanted <- c(wanted, list(list(db = db, id = id)))
  }
  # Plus the top curated sets by pooled p-value, then |NES| — "top-enriching" read off
  # the pooled correction, so the choice is not made under a laxer test than the one
  # the tables report. Selecting from `sweep_df` also means an alias copy can never be
  # chosen, so a re-pinned Hallmark set gets one substrate file under its canonical
  # collection instead of two identical curves under two database names.
  # The quota is split EVENLY BETWEEN DIRECTIONS: ranked on pooled
  # p-value alone, all five slots went to translation and ribosome sets on the
  # paired-blood side, which leaves the downstream running-sum comparison with no
  # up-going curated comparator for the mouse arms — the one thing it is built to do.
  top <- sweep_df |>
    dplyr::filter(population == pop, database %in% CURATED_DBS) |>
    dplyr::group_by(direction) |>
    dplyr::arrange(padj_pooled, pvalue, dplyr::desc(abs(nes)), .by_group = TRUE) |>
    dplyr::slice_head(n = RUNSUM_HALF) |>
    dplyr::ungroup()
  wanted <- c(wanted, lapply(seq_len(nrow(top)),
                             function(i) list(db = top$database[i], id = top$pathway_id[i])))
  # A named set can also rank into the top-N quota, and the write would then be issued
  # twice: harmless on disk (the second overwrites the first) but it would put two rows
  # for one curve into the index, and the downstream figure loop reads the index.
  wanted <- wanted[!duplicated(vapply(wanted, function(w) paste(w$db, w$id), character(1)))]

  for (w in wanted) {
    g <- gsea_cells[[paste(pop, w$db)]]
    genes <- COLLECTIONS[[w$db]]$sets[[w$id]]
    if (is.null(genes)) next
    if (!w$id %in% g@result$ID) {
      message(sprintf("    %s x %s: not scored (effective size outside [%d,%d]) — no substrate",
                      pop, w$id, MINSZ, MAXSZ))
      next
    }
    df <- runsum_table(g, w$id, genes, pop)
    # Set ids come from curated collections and CollecTRI TF names; anything outside
    # a filename-safe alphabet is substituted so the path can never be ambiguous. The
    # untouched id travels in the `gene_set` column and in the index, so the mapping
    # from file to set is always recoverable.
    fp <- file.path(TBL, sprintf("runsum_interactive_%s_%s.csv", POPS[[pop]],
                                 gsub("[^A-Za-z0-9._-]", "_", w$id)))
    utils::write.csv(df, fp, row.names = FALSE)
    r <- g@result[g@result$ID == w$id, ]
    runsum_index[[length(runsum_index) + 1L]] <- tibble::tibble(
      file = basename(fp), population = pop, database = w$db, gene_set = w$id,
      contrast = contrast_label(pop), n_ranked = nrow(df), set_size = r$setSize[1],
      nes = r$NES[1], padj_in_database = r$p.adjust[1],
      padj_pooled = sweep_df$padj_pooled[sweep_df$population == pop & sweep_df$pathway_id == w$id &
                                        sweep_df$database == w$db][1],
      # TRUE when this curve was emitted because the set was NAMED — every mouse-derived
      # arm, plus the config's runsum_always list — rather than because it ranked into the
      # top-N quota. A reader of the index can then tell a guaranteed comparator from one
      # that happens to be top-ranked in this population and may vanish in the next.
      always_emitted = identical(w$db, GATE_DB) || w$id %in% RUNSUM_ALWAYS)
  }
}
runsum_idx <- dplyr::bind_rows(runsum_index)
if (nrow(runsum_idx) == 0)
  stop("[14] no running-sum substrate written at all — the downstream comparison has no input.")
readr::write_csv(round_numeric_cols(runsum_idx), file.path(TBL, "runsum_interactive_index.csv"))
message(sprintf("  %d running-sum substrate table(s) written", nrow(runsum_idx)))

# Schema guard: the downstream widget joins on these names in this order.
RUNSUM_COLS <- c("rank", "gene", "stat", "running_es", "hit", "leading_edge",
                 "gene_set", "population", "contrast")
ref_runsum <- file.path(RESULTS, "05_scoring", CFG$paths$stage_tables_subdir %||% "tables",
                        "runsum_interactive_treg_WT_heat_up.csv")
if (file.exists(ref_runsum)) {
  ref_cols <- colnames(readr::read_csv(ref_runsum, n_max = 1, show_col_types = FALSE,
                                       progress = FALSE))
  if (!identical(ref_cols, RUNSUM_COLS))
    stop("[14] the published stage-05 runsum schema is ", paste(ref_cols, collapse = ","),
         " but this stage writes ", paste(RUNSUM_COLS, collapse = ","),
         " — the downstream widget joins on these, so fix before proceeding.")
  mine <- colnames(readr::read_csv(file.path(TBL, runsum_idx$file[1]), n_max = 1,
                                   show_col_types = FALSE, progress = FALSE))
  stopifnot("emitted runsum schema differs from the published one" = identical(mine, ref_cols))
  message("  runsum schema matches the published stage-05 tables exactly.")
}

# ============================================================================
# 9. PROGENy — decoupleR MLM, no gene-set list involved
# ============================================================================
## Method-orthogonal corroboration. PROGENy footprints are continuous-weight models
## of fourteen signalling pathways fitted to perturbation experiments, so this arm
## does not inherit the size floors, curation choices or set-overlap structure the
## fgsea sweep does. It is run two ways:
##   (a) on the moderated-t CONTRAST statistics, mirroring the mouse anchor exactly —
##       one column per population, so each score is the SF-versus-PB shift;
##   (b) on the donor-level pseudobulk, which the contrast statistics cannot give:
##       a per-donor activity that can be tested paired, so the (a) score is not the
##       only evidence for a direction.
## `.mor = "weight"` because PROGENy weights are continuous. Passing "mor" would
## silently binarise them.

message("[9] PROGENy pathway activity (decoupleR MLM, human model) ...")
if (!requireNamespace("decoupleR", quietly = TRUE))
  stop("[14] decoupleR is required. BiocManager::install('decoupleR').")
if (!requireNamespace("progeny", quietly = TRUE))
  stop("[14] progeny is required. BiocManager::install('progeny').")

PG      <- UE$progeny
PG_ORG  <- PG$organism %||% "Human"
PG_TOP  <- as.integer(PG$top %||% 500L)
PG_MIN  <- as.integer(PG$minsize %||% 5L)

pg_mat <- progeny::getModel(PG_ORG, top = PG_TOP)
net_progeny <- data.frame(
  target = rep(rownames(pg_mat), times = ncol(pg_mat)),
  source = rep(colnames(pg_mat), each = nrow(pg_mat)),
  weight = as.numeric(as.matrix(pg_mat)),
  stringsAsFactors = FALSE)
net_progeny <- net_progeny[net_progeny$weight != 0, ]
stopifnot("PROGENy model must carry 14 pathways" = length(unique(net_progeny$source)) == 14L)
message(sprintf("  PROGENy %s top=%d: %d pathways, %d weighted edges, %d target genes",
                PG_ORG, PG_TOP, length(unique(net_progeny$source)), nrow(net_progeny),
                length(unique(net_progeny$target))))

# --- (a) contrast statistics -------------------------------------------------
## The three populations do not share a gene universe: filterByExpr kept a slightly
## different set for each. A gene absent from one population's DE table is padded
## with t = 0 ("no evidence of a shift"), which is what the mouse anchor does with
## NA, and the padding count is published per population so the padding is visible
## rather than assumed harmless.
de_stats <- lapply(names(POPS), function(pop) {
  p <- file.path(RANKED_DIR, sprintf("de_SFvsPB_%s.csv", POPS[[pop]]))
  if (!file.exists(p)) stop("[14] DE table not found: ", p)
  d <- readr::read_csv(p, show_col_types = FALSE, progress = FALSE)
  d <- d[!is.na(d$gene_symbol) & nzchar(d$gene_symbol) & !is.na(d$stat), ]
  d <- d[order(-abs(d$stat)), ]
  d <- d[!duplicated(d$gene_symbol), ]
  stats::setNames(d$stat, d$gene_symbol)
})
names(de_stats) <- names(POPS)
uni <- sort(unique(unlist(lapply(de_stats, names), use.names = FALSE)))
tmat <- vapply(de_stats, function(v) { out <- rep(0, length(uni)); names(out) <- uni
                                       out[names(v)] <- v; out }, numeric(length(uni)))
rownames(tmat) <- uni
pad <- vapply(de_stats, function(v) length(uni) - length(v), integer(1))
message(sprintf("  contrast t-matrix: %d genes x %d populations (zero-padded: %s)",
                nrow(tmat), ncol(tmat),
                paste(sprintf("%s=%d", names(pad), pad), collapse = " ")))

progeny_contrast_raw <- decoupleR::run_mlm(
  mat = tmat, network = net_progeny, .source = "source", .target = "target",
  .mor = "weight", minsize = PG_MIN)

footprint <- net_progeny |>
  dplyr::filter(target %in% rownames(tmat)) |>
  dplyr::group_by(source) |>
  dplyr::summarise(set_size = dplyr::n(),
                   core_enrichment = paste(sort(unique(target)), collapse = "/"),
                   .groups = "drop")

## NAMING, FLAGGED FOR WHOEVER REGENERATES THIS STAGE. The activity statistic is
## written into a column called `nes` for schema compatibility with
## master_gsea_table and with the ten fgsea outputs it lands beside — but it is NOT a
## normalized enrichment score and shares no scale with them. These 42 rows span
## -3.87 to +10.85 where the 33,955-row fgsea sweep spans -3.76 to +2.75, so 15 of
## them exceed the largest absolute NES anywhere in the sweep; JAK-STAT at +9.4 in
## Treg is an ordinary MLM score and an impossible NES. The stage README carries the
## warning where a reader of the CSV will meet it. If the schema constraint is ever
## relaxed, rename this to `mlm_score` (or add a `statistic_kind` column) and update
## that caption in the same commit — the ambiguity is deliberate, not a defect.
progeny_activity <- progeny_contrast_raw |>
  dplyr::filter(!is.na(score)) |>
  dplyr::group_by(condition) |>
  dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::transmute(pathway_id = paste0("PROGENY_", source), pathway_name = source,
                   database = "PROGENy", population = condition,
                   contrast = vapply(condition, contrast_label, character(1)),
                   nes = score, pvalue = p_value, padj = padj,
                   direction = ifelse(score > 0, "up", "down")) |>
  dplyr::left_join(dplyr::rename(footprint, pathway_name = source), by = "pathway_name") |>
  dplyr::select(pathway_id, pathway_name, database, population, contrast, nes, pvalue,
                padj, set_size, core_enrichment, direction) |>
  dplyr::arrange(population, pvalue)
readr::write_csv(round_numeric_cols(progeny_activity), file.path(TBL, "progeny_activity.csv"))
message(sprintf("  progeny_activity.csv: %d rows (%d pathways x %d populations)",
                nrow(progeny_activity), length(unique(progeny_activity$pathway_name)),
                length(unique(progeny_activity$population))))

# --- (b) donor-level pseudobulk ---------------------------------------------
## The counts are Ensembl-keyed and gene_symbols.csv exists precisely for this join.
## Duplicated symbols are SUMMED (they are counts from distinct Ensembl ids of one
## gene, so summing is the only aggregation that keeps a count a count — the max-|t|
## collapse used on the ranked lists is a statistic rule and does not apply here).
## The retained gene universe is that population's own DE-table universe, so no new
## expression filter is introduced here and the two PROGENy arms see the same genes.
## log2 CPM is row-centred WITHIN each population, so a score compares donors of one
## sorted population and never carries a between-population difference.

counts_path  <- file.path(RANKED_DIR, "pseudobulk_counts.csv")
coldata_path <- file.path(RANKED_DIR, "pseudobulk_coldata.csv")
genemap_path <- file.path(RANKED_DIR, "gene_symbols.csv")
stopifnot("pseudobulk_counts.csv missing"  = file.exists(counts_path),
          "pseudobulk_coldata.csv missing" = file.exists(coldata_path),
          "gene_symbols.csv missing — the Python->R seam map is mandatory" =
            file.exists(genemap_path))

counts_df <- readr::read_csv(counts_path, show_col_types = FALSE, progress = FALSE)
coldata   <- readr::read_csv(coldata_path, show_col_types = FALSE, progress = FALSE)
genemap   <- readr::read_csv(genemap_path, show_col_types = FALSE, progress = FALSE)
colnames(counts_df)[1] <- "sample_id"
colnames(coldata)[1]   <- "sample_id"

cmat <- t(as.matrix(counts_df[, -1, drop = FALSE]))
colnames(cmat) <- counts_df$sample_id
sym <- genemap$gene_symbol[match(rownames(cmat), genemap$ensembl_id)]
keep <- !is.na(sym) & nzchar(sym)
message(sprintf("  donor pseudobulk: %d x %d; %d of %d Ensembl ids carry a symbol",
                nrow(cmat), ncol(cmat), sum(keep), nrow(cmat)))
if (mean(keep) < 0.5)
  stop("[14] fewer than half the pseudobulk Ensembl ids map to a symbol — check gene_symbols.csv.")
cmat <- cmat[keep, , drop = FALSE]
cmat <- rowsum(cmat, group = sym[keep])

TISSUE_NUM <- CFG$design$tissue_levels$synovial_fluid   %||% "synovial_fluid"
TISSUE_DEN <- CFG$design$tissue_levels$peripheral_blood %||% "peripheral_blood"

## These two tables escape the schema constraint noted above and so name their
## statistics honestly: `activity` here and `mean_difference` in the paired table, both
## in decoupleR MLM units. Neither is an enrichment score, and neither shares a scale
## with any `gsea_*.csv` column.
donor_rows <- list(); paired_rows <- list()
for (pop in names(POPS)) {
  cd <- coldata[coldata$coarse_label == pop, ]
  missing_cols <- setdiff(cd$sample_id, colnames(cmat))
  if (length(missing_cols))
    stop("[14] coldata names pseudobulk samples absent from the count matrix: ",
         paste(missing_cols, collapse = ", "))
  m  <- cmat[intersect(rownames(cmat), names(de_stats[[pop]])), cd$sample_id, drop = FALSE]
  lib <- colSums(m)
  if (any(lib == 0))
    stop("[14] a pseudobulk sample has zero counts over this population's DE universe: ",
         paste(colnames(m)[lib == 0], collapse = ", "))
  lg  <- log2(t(t(m) / lib) * 1e6 + 1)
  lg  <- lg - rowMeans(lg)                      # centre within this population
  # A gene that is constant across this population's samples carries no information
  # after centring and would enter the model as an all-zero column.
  lg  <- lg[apply(lg, 1, function(x) any(x != 0)), , drop = FALSE]
  act <- decoupleR::run_mlm(mat = lg, network = net_progeny, .source = "source",
                            .target = "target", .mor = "weight", minsize = PG_MIN)
  act <- dplyr::left_join(act, dplyr::select(cd, sample_id, donor, tissue, coarse_label,
                                             n_cells),
                          by = c("condition" = "sample_id"))
  donor_rows[[pop]] <- act |>
    dplyr::transmute(population = pop, sample_id = condition, donor = donor, tissue = tissue,
                     n_cells = n_cells, pathway_name = source, activity = score,
                     pvalue = p_value, n_genes_scored = nrow(lg))

  # Paired SF-vs-PB per pathway on donors present in BOTH arms. Paired because the
  # design is paired; an unpaired test here would throw away the blocking the DE
  # model already used and is not the same question.
  d <- donor_rows[[pop]]
  both <- intersect(d$donor[d$tissue == TISSUE_NUM], d$donor[d$tissue == TISSUE_DEN])
  key <- paste(d$pathway_name, d$donor, d$tissue)
  for (pw in unique(d$pathway_name)) {
    # Index BY DONOR so the two vectors are aligned donor-for-donor. Filtering each
    # arm separately and trusting row order would silently pair the wrong donors when
    # an arm is missing one, which is exactly the case here (six donors span both arms
    # in Treg, seven in Tcon and CD8).
    sf <- d$activity[match(paste(pw, both, TISSUE_NUM), key)]
    pb <- d$activity[match(paste(pw, both, TISSUE_DEN), key)]
    ok <- !is.na(sf) & !is.na(pb)
    tt <- if (sum(ok) >= 3) stats::t.test(sf[ok], pb[ok], paired = TRUE) else NULL
    paired_rows[[paste(pop, pw)]] <- tibble::tibble(
      population = pop, contrast = contrast_label(pop), pathway_name = pw,
      n_paired_donors = sum(ok),
      mean_sf = mean(sf[ok]), mean_pb = mean(pb[ok]),
      mean_difference = mean(sf[ok] - pb[ok]),
      t_statistic = if (is.null(tt)) NA_real_ else unname(tt$statistic),
      pvalue = if (is.null(tt)) NA_real_ else tt$p.value,
      ci_low = if (is.null(tt)) NA_real_ else tt$conf.int[1],
      ci_high = if (is.null(tt)) NA_real_ else tt$conf.int[2])
  }
}

progeny_donor <- dplyr::bind_rows(donor_rows) |> dplyr::arrange(population, pathway_name, donor)
readr::write_csv(round_numeric_cols(progeny_donor), file.path(TBL, "progeny_donor_activity.csv"))

progeny_paired <- dplyr::bind_rows(paired_rows) |>
  dplyr::group_by(population) |>
  dplyr::mutate(padj = stats::p.adjust(pvalue, method = "BH"),
                direction = ifelse(is.na(mean_difference), NA_character_,
                                   ifelse(mean_difference > 0, "up", "down"))) |>
  dplyr::ungroup() |>
  dplyr::arrange(population, pvalue)
readr::write_csv(round_numeric_cols(progeny_paired), file.path(TBL, "progeny_sf_vs_pb.csv"))

saveRDS(list(contrast = progeny_contrast_raw, net = net_progeny,
             donor = progeny_donor, paired = progeny_paired),
        file.path(DIR_OBJ, "14_progeny.rds"))
message(sprintf("  progeny_donor_activity.csv: %d rows; progeny_sf_vs_pb.csv: %d rows",
                nrow(progeny_donor), nrow(progeny_paired)))

# ============================================================================
# 10. ACCEPTANCE CHECKS + headline log
# ============================================================================

message("[10] acceptance checks ...")
must_exist <- c("geneset_manifest.csv", "ranked_list_keycheck.csv", "geneset_overlap.csv",
                "geneset_alias_map.csv", "wt_heat_up_reproduction.csv", "gsea_all.csv",
                "gsea_pooled_summary_by_db.csv", "runsum_interactive_index.csv",
                "progeny_activity.csv", "progeny_donor_activity.csv", "progeny_sf_vs_pb.csv")
for (f in must_exist)
  stopifnot(setNames(file.exists(file.path(TBL, f)), paste(f, "missing")))
for (pop in names(POPS)) for (db in DB_ORDER)
  stopifnot(setNames(file.exists(file.path(TBL, sprintf("gsea_%s_%s.csv", POPS[[pop]], db))),
                     sprintf("gsea_%s_%s.csv missing", POPS[[pop]], db)))

## THE REGRESSION GUARD FOR THIS STAGE'S PAST DEFECT. The pooled family must contain
## each hypothesis once. When it did not, the six re-pinned Hallmark sets entered BH
## twice per population and the rank denominators quoted downstream (11,242 / 11,465 /
## 11,248) were each six too high — a wrong number that no error surfaced. Asserted, not
## logged, because a duplicate here silently corrupts every pooled q-value and every
## "rank k of N" a reader takes from this table.
dup_check <- sweep_out |>
  dplyr::count(population, pathway_id) |>
  dplyr::filter(n > 1)
if (nrow(dup_check) > 0)
  stop(sprintf(paste0("[14] gsea_all.csv carries %d duplicated (population, pathway_id) pair(s), ",
                      "e.g. %s in %s. The pooled BH family would contain the same hypothesis ",
                      "more than once and every rank denominator would be inflated. Section 2b's ",
                      "alias resolution did not cover this collision."),
               nrow(dup_check), dup_check$pathway_id[1], dup_check$population[1]))
message(sprintf("  gsea_all.csv holds one row per (population, pathway_id): %s",
                paste(sprintf("%s %d", names(POPS),
                              vapply(names(POPS),
                                     function(p) sum(sweep_out$population == p), integer(1))),
                      collapse = ", ")))
## And the alias copies must agree with the kept copies to fgsea's permutation noise. A
## larger gap would mean the two copies are not the same test, so dropping one would be
## discarding evidence rather than removing a duplicate.
if (nrow(ALIAS) > 0) {
  worst <- max(alias_out$abs_nes_difference, na.rm = TRUE)
  if (worst > NES_TOL)
    stop(sprintf(paste0("[14] a dropped alias copy differs from the copy kept in the pooled ",
                        "family by |NES| %.4f, above the %g tolerance fgsea's stochastic ",
                        "normalisation explains. The two copies are then not the same test and ",
                        "one of them is not a duplicate — do not deduplicate until this is ",
                        "understood."), worst, NES_TOL))
  message(sprintf("  alias copies agree with the kept copies to |NES| <= %.2e (tolerance %g).",
                  worst, NES_TOL))
}
## Pooling over a wider family almost always LOOSENS a row's significance, but not as
## a theorem: BH divides by rank as well as multiplying by family size, so a row whose
## rank rises faster than the family grows can come out marginally tighter. The count
## is reported rather than asserted, so a surprising number is visible instead of
## either crashing the run or passing unseen.
n_tighter <- sum(sweep_df$padj_pooled < sweep_df$padj - 1e-12, na.rm = TRUE)
message(sprintf("  %d of %d rows have padj_pooled below their per-database padj (%.3f%%)",
                n_tighter, nrow(sweep_df), 100 * n_tighter / nrow(sweep_df)))

message("")
message("  --- top 15 pooled-significant sets, Treg ---")
top_treg <- sweep_out |>
  dplyr::filter(population == "Treg", padj_pooled < FDR) |>
  dplyr::arrange(padj_pooled, dplyr::desc(abs(nes))) |>
  dplyr::slice_head(n = 15)
for (i in seq_len(nrow(top_treg))) {
  r <- top_treg[i, ]
  message(sprintf("  %2d. %-14s %-58s NES %+.2f  padj_pooled %.2e  n=%d",
                  i, r$database, substr(r$pathway_id, 1, 58), r$nes, r$padj_pooled, r$set_size))
}
message("")
message("  --- where the mouse-derived up arms rank inside their own population's sweep ---")
for (pop in names(POPS)) {
  s <- sweep_out[sweep_out$population == pop, ]
  s <- s[order(s$padj_pooled, -abs(s$nes)), ]
  for (sig in names(COLLECTIONS[[GATE_DB]]$sets)) {
    k <- which(s$pathway_id == sig & s$database == GATE_DB)
    if (!length(k)) { message(sprintf("  %-4s %-16s not scored", pop, sig)); next }
    message(sprintf("  %-4s %-16s NES %+.2f  padj_pooled %.2e  rank %d of %d;  %d sets carry a larger |NES| at FDR<%.2g",
                    pop, sig, s$nes[k], s$padj_pooled[k], k, nrow(s),
                    sum(abs(s$nes) > abs(s$nes[k]) & s$padj_pooled < FDR, na.rm = TRUE), FDR))
  }
}
message("")
message("  --- PROGENy on the contrast statistics ---")
for (i in seq_len(nrow(progeny_activity))) {
  r <- progeny_activity[i, ]
  if (r$padj < FDR)
    message(sprintf("  %-4s %-10s score %+.2f  padj %.2e  (%s)", r$population, r$pathway_name,
                    r$nes, r$padj, r$direction))
}

message("=================================================================")
message("14_unbiased_enrichment COMPLETE")
message("  Tables:  ", TBL)
message("  Objects: ", file.path(DIR_OBJ, "14_genesets.rds"), ", ", DIR_GSEA, "/, ",
        file.path(DIR_OBJ, "14_progeny.rds"))
message("  Run 14_unbiased_enrichment_viz.R for figures.")
message("=================================================================")
