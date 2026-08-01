#!/usr/bin/env Rscript
# =============================================================================
# 15_coresh_search.R  --  COMPUTE (no plots)
# =============================================================================
# Co-regulation search of the public HUMAN GEO compendium for the JIA synovial-
# fluid-versus-paired-blood niche contrast, and pre-ranked GSEA of the modules it
# returns on the same ranked list they were seeded from.
#
# The question. The mouse anchor asked CoReSh which public co-regulation
# neighbourhoods its 39 degC contrast sits in. Nothing had asked the same of the
# human niche contrast. This script asks it: across ~44,000 public human GEO
# datasets, in which ones do the genes that rise in the inflamed synovial niche
# co-vary, and what else moves with them there?
#
# What a "module" is here. CoReSh scores each public dataset by how much of its
# variance the query genes jointly explain (pctVar, a PCA-inspired quantity). For
# a top-ranked dataset, projecting every gene onto the query direction gives
# gene-level loadings; the top-|loading| genes are that dataset's co-regulation
# partners of the query. The resulting set is a DATA-DRIVEN CO-REGULATION
# NEIGHBOURHOOD mined from public variance structure -- not a curated ontology
# term, and not a claim about mechanism. It is named for how it was derived:
# CORESH_<population>_up_<gate>_<GSE>.
#
# Circularity, stated up front. A module seeded from a query contains that
# query's genes by construction, so its enrichment on the seeding ranked list is
# partly guaranteed. The `seeded_from_this_population` column marks exactly which
# rows carry that circularity, and the fraction of each module that is seed
# rather than newly recruited gene is published as `frac_seed_genes`. Read the
# enrichment as a description of what public biology the niche signature co-moves
# with, never as independent evidence for it.
#
# Query construction MIRRORS the mouse anchor (mouse_anchor/02_analysis/scripts/
# 07_coresh_search.R + its coresh.query_signatures block): the UP arm of the
# contrast at two stringency gates, fdr_only and fdr_logfc. Here the contrast is
# the JIA SF-vs-PB contrast and the gates are re-derived from the same frozen
# limma-voom DE table that produced ranked_<population>.tsv.
#
# The ID seam. Ranked lists and every reference gene set are keyed by HGNC
# SYMBOL; the compendium chunks are keyed by INTEGER ENTREZ ID. match() on a
# character vector against an integer vector returns all-NA silently, so a
# symbol-keyed query would return size ~= 0 and pctVar noise with no error. The
# script maps symbol -> Entrez explicitly, publishes the mapping loss, and
# refuses to run if the realised per-dataset overlap collapses.
#
# Inputs (read-only):
#   03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv     signed-t ranked lists
#   03_results/03_pseudobulk/tables/de_SFvsPB_{treg,tcon,cd8}.csv  the gates come from here
#   $CORESH_CHUNKS/hsa/*_full_objects.qs2                          shared read-only cache
#
# Outputs (03_results/15_coresh_search/tables/):
#   coresh_query_provenance.csv          what was queried + Entrez mapping loss
#   coresh_hits.csv                      the ranked compendium, one row per dataset per query
#   coresh_derived_sets.csv              the derived modules, with sizes and members
#   coresh_derived_sets.gmt              the same modules in GMT form
#   coresh_derived_gsea.csv              fgsea of those modules on the ranked lists
#   coresh_derived_annotation.csv        what each recovered dataset is, from the cache itself
#   coresh_module_sizes.csv              module-size / overlap distribution
#   runsum_interactive_coresh_<pop>_<set>.csv   running-sum substrate for interactive curves
# Checkpoints (03_results/objects/): coresh_hsa_*.rds
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/15_coresh_search.R
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(qs2)
  library(clusterProfiler)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

YAML_CONFIG <- yaml::read_yaml("02_analysis/config/analysis_config.yaml")

STAGE       <- "15_coresh_search"
SCRIPT_PATH <- "02_analysis/scripts/15_coresh_search.R"
RESULTS     <- YAML_CONFIG$paths$results %||% "03_results/"
DIR_OBJECTS <- YAML_CONFIG$paths$objects %||% "03_results/objects/"
TBL_DIR     <- file.path(RESULTS, STAGE, "tables")
OVERVIEW    <- YAML_CONFIG$figures$overview_dir %||% "_overview"
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TBL_DIR, OVERVIEW), recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)

## Declared stage guard -- the stage must exist in analysis_config.yaml:stages.
stage_ids <- vapply(YAML_CONFIG$stages, function(s) s$id %||% "", character(1))
if (!STAGE %in% stage_ids)
  stop("[15_coresh] stage '", STAGE, "' is not declared in analysis_config.yaml:stages.")

CFG <- YAML_CONFIG$coresh %||% list()
if (length(CFG) == 0L)
  stop("[15_coresh] no `coresh:` block in analysis_config.yaml -- nothing is configured.")

SPECIES        <- CFG$species             %||% "human"
CHUNK_SUBDIR   <- CFG$chunks_subdir       %||% "hsa"
TOP_N          <- as.integer(CFG$top_n_hits      %||% 5L)
MIN_Q          <- as.integer(CFG$min_query_size  %||% 3L)
N_CORES        <- as.integer(CFG$n_cores         %||% 4L)
USE_PVALUES    <- isTRUE(CFG$pvalues)
N_DERIVE       <- as.integer(CFG$n_derive        %||% 50L)
JACCARD        <- as.numeric(CFG$jaccard         %||% 0.8)
MIN_SET        <- as.integer(CFG$derived_min_size %||% 15L)
MAX_SET        <- as.integer(CFG$derived_max_size %||% 500L)
N_TERMS        <- as.integer(CFG$n_annotation_terms %||% 15L)
PRIMARY_POP    <- CFG$primary_population  %||% "treg"

THR            <- YAML_CONFIG$thresholds %||% list()
DE_FDR         <- as.numeric(THR$de_fdr        %||% 0.05)
DE_LOGFC       <- as.numeric(THR$de_logfc      %||% 1.0)
GSEA_MIN       <- as.integer(THR$gsea_min_size %||% 5L)
GSEA_MAX       <- as.integer(THR$gsea_max_size %||% 500L)
GSEA_SEED      <- as.integer(THR$gsea_seed     %||% 123L)
GSEA_NPERM     <- as.integer(THR$gsea_nperm    %||% 100000L)
GSEA_FDR       <- as.numeric(THR$gsea_fdr      %||% 0.05)
RUNSUM_TOP     <- as.integer(YAML_CONFIG$figures$running_sum_top %||% 5L)

## hsa chunks demand HUMAN Entrez. A mouse query against hsa (or vice versa)
## returns near-zero sizes with NO error -- refuse rather than silently mis-run.
if (!identical(SPECIES, "human"))
  stop("[15_coresh] coresh.species must be 'human' for the hsa compendium; got '", SPECIES, "'.")
if (!identical(CHUNK_SUBDIR, "hsa"))
  stop("[15_coresh] coresh.chunks_subdir must be 'hsa' for human data; got '", CHUNK_SUBDIR, "'.")

message(sprintf(
  "[15_coresh] species=%s subdir=%s top_n_hits=%d n_cores=%d pvalues=%s | derived: n=%d size=[%d,%d] jaccard=%.2f",
  SPECIES, CHUNK_SUBDIR, TOP_N, N_CORES, USE_PVALUES, N_DERIVE, MIN_SET, MAX_SET, JACCARD))

# =============================================================================
# 1. Canonical CoReSh skill scripts (sym2ent FIRST -- the bridge asserts on it)
# =============================================================================

CORESH_LIB <- "01_modules/SciAgent-toolkit/skills/coresh-signature-search/scripts"
if (!dir.exists(CORESH_LIB))
  stop("[15_coresh] CoReSh skill scripts not found at ", CORESH_LIB,
       " -- the coresh-signature-search skill must be present under 01_modules/SciAgent-toolkit.")
source(file.path(CORESH_LIB, "symbols_to_entrez.R"))      # sym2ent() / ent2sym()
source(file.path(CORESH_LIB, "coresh_batch.R"))           # coresh_batch() / coreshMatch()
source(file.path(CORESH_LIB, "extract_gene_loadings.R"))  # build_coresh_gmt() -- needs sym2ent

# =============================================================================
# 2. Pre-flight: the shared read-only chunk cache (STOP loudly; never fabricate)
# =============================================================================
# Resolution order: $CORESH_CHUNKS (the container convention) > coresh.chunks_fallback.
# The species subdir is appended if the resolved dir does not already end in it.
# The analysis container never sees a Synapse token -- it consumes pre-cached
# chunks only, so an absent cache is a provisioning problem, not an analysis one.

resolve_chunk_dir <- function() {
  base <- Sys.getenv("CORESH_CHUNKS", unset = "")
  if (!nzchar(base)) base <- CFG$chunks_fallback %||%
    "00_data/references/coresh/current/preprocessed_chunks"
  cand <- if (basename(base) == CHUNK_SUBDIR) base else file.path(base, CHUNK_SUBDIR)
  if (!dir.exists(cand))
    stop("[15_coresh] CoReSh ", CHUNK_SUBDIR, " chunk directory not found: ", cand,
         "\n  The compendium is staged out-of-band into the shared reference cache",
         "\n  (Synapse syn66227307, ~20 GB); it cannot be downloaded from this container.",
         "\n  Set $CORESH_CHUNKS to the preprocessed_chunks dir and re-run.")
  normalizePath(cand, mustWork = TRUE)
}
CHUNK_DIR   <- resolve_chunk_dir()
CHUNK_FILES <- list.files(CHUNK_DIR, pattern = "_full_objects\\.qs2$", full.names = TRUE)
if (length(CHUNK_FILES) == 0L)
  stop("[15_coresh] no *_full_objects.qs2 chunks under ", CHUNK_DIR,
       " -- the directory exists but is empty or mis-shaped (expect ~89 hsa chunks).")

## Cache identity, so the ranking is attributable to a specific snapshot.
manifest_fp <- file.path(dirname(dirname(CHUNK_DIR)), "MANIFEST.json")
CACHE_TAG   <- NA_character_; CACHE_DATE <- NA_character_; CACHE_SYN <- NA_character_
if (file.exists(manifest_fp) && requireNamespace("jsonlite", quietly = TRUE)) {
  mf <- jsonlite::fromJSON(manifest_fp)
  CACHE_TAG  <- as.character(mf$snapshot_tag  %||% NA)
  CACHE_DATE <- as.character(mf$downloaded_at %||% NA)
  CACHE_SYN  <- as.character(mf$synapse_id    %||% NA)
}
message(sprintf("[15_coresh] compendium: %d %s chunk(s) at %s (snapshot %s, %s)",
                length(CHUNK_FILES), CHUNK_SUBDIR, CHUNK_DIR, CACHE_TAG, CACHE_DATE))

# =============================================================================
# 3. Build one query per configured (population, gate) entry
# =============================================================================
# The gate is applied to the frozen limma-voom DE table; the ranked list is the
# same fit collapsed to one row per symbol, so every gated symbol is guaranteed
# to be present in the ranked list. That is asserted, not assumed.

QSPEC <- CFG$query_signatures %||% list()
if (is.null(QSPEC$sets) || length(QSPEC$sets) == 0L)
  stop("[15_coresh] coresh.query_signatures.sets must be a non-empty list.")
DE_TMPL  <- QSPEC$de_table %||% "03_results/03_pseudobulk/tables/de_SFvsPB_{population}.csv"
RNK_TMPL <- QSPEC$ranked   %||% "03_results/03_pseudobulk/tables/ranked_{population}.tsv"
fill_tmpl <- function(tmpl, population) gsub("{population}", population, tmpl, fixed = TRUE)

POPULATIONS <- unique(vapply(QSPEC$sets, function(s) as.character(s$population), character(1)))
CONTRAST_LABEL <- c(treg = "SF_vs_PB_Treg", tcon = "SF_vs_PB_Tcon", cd8 = "SF_vs_PB_CD8")

## --- ranked lists (HGNC symbol -> signed moderated t, descending) -----------
read_ranked <- function(population) {
  fp <- fill_tmpl(RNK_TMPL, population)
  if (!file.exists(fp)) stop("[15_coresh] ranked list not found: ", fp)
  r <- utils::read.table(fp, sep = "\t", header = FALSE,
                         col.names = c("symbol", "stat"), quote = "", comment.char = "")
  r <- r[!is.na(r$stat) & nzchar(r$symbol), ]
  r <- r[!duplicated(r$symbol), ]
  ## The documented silent failure: an Ensembl-keyed list intersects every
  ## symbol-keyed reference at ~zero and fgsea reports it as empty, not as error.
  ens_frac <- mean(grepl("^ENSG[0-9]+", r$symbol))
  if (ens_frac > 0.5)
    stop("[15_coresh] ", fp, " looks Ensembl-keyed (", round(100 * ens_frac), "% ENSG ids). ",
         "Every reference set here matches on HGNC symbol; re-key the ranked list and re-run.")
  v <- sort(stats::setNames(r$stat, r$symbol), decreasing = TRUE)
  v
}
RANKED <- lapply(stats::setNames(POPULATIONS, POPULATIONS), read_ranked)
for (p in POPULATIONS)
  message(sprintf("  ranked_%s: %d symbols, t in [%.2f, %.2f]",
                  p, length(RANKED[[p]]), min(RANKED[[p]]), max(RANKED[[p]])))

## --- gate the DE table, then map symbol -> Entrez, accounting for every loss --
build_query <- function(spec) {
  population <- as.character(spec$population)
  gate       <- as.character(spec$gate)
  direction  <- as.character(spec$direction %||% "up")
  if (!identical(direction, "up"))
    stop("[15_coresh] only the up arm is queried (mirrors the mouse anchor); got '", direction, "'.")

  de_fp <- fill_tmpl(DE_TMPL, population)
  if (!file.exists(de_fp)) stop("[15_coresh] DE table not found: ", de_fp)
  de <- utils::read.csv(de_fp)
  need <- c("gene_symbol", "log2FoldChange", "padj")
  if (!all(need %in% names(de)))
    stop("[15_coresh] ", de_fp, " missing column(s): ",
         paste(setdiff(need, names(de)), collapse = ", "))

  keep <- !is.na(de$padj) & de$padj < DE_FDR & de$log2FoldChange > 0
  if (identical(gate, "fdr_logfc")) keep <- keep & de$log2FoldChange >= DE_LOGFC
  else if (!identical(gate, "fdr_only"))
    stop("[15_coresh] unknown gate '", gate, "' (expected fdr_only or fdr_logfc).")
  symbols <- unique(de$gene_symbol[keep])
  symbols <- symbols[nzchar(symbols) & !is.na(symbols)]
  if (length(symbols) == 0L)
    stop("[15_coresh] gate ", gate, " on ", de_fp, " selected no genes.")

  ## The gate is taken from the same fit the ranked list came from, so the query
  ## must be a subset of the ranked list. If it is not, the two have drifted.
  not_ranked <- setdiff(symbols, names(RANKED[[population]]))
  if (length(not_ranked) > 0L)
    stop("[15_coresh] ", length(not_ranked), " gated symbol(s) absent from ranked_", population,
         ".tsv (e.g. ", paste(utils::head(not_ranked, 5), collapse = ", "),
         ") -- DE table and ranked list have drifted apart.")

  ## symbol -> Entrez. mapIds is called directly so the LOSS can be counted; the
  ## result is then asserted identical to the skill helper's own output.
  mapped <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
    column = "ENTREZID", multiVals = "first"))
  unmapped   <- symbols[is.na(mapped)]
  entrez_all <- as.integer(stats::na.omit(mapped))
  entrez     <- unique(entrez_all)          # a duplicated id would inflate k and double-count
  via_skill  <- suppressWarnings(sym2ent(symbols, species = SPECIES))
  stopifnot(identical(sort(unique(via_skill)), sort(entrez)))

  qname <- sprintf("%s_up_%s", population, gate)
  message(sprintf("  %-20s gate=%-9s %5d symbols -> %5d Entrez (%d unmapped, %d duplicate ids collapsed)",
                  qname, gate, length(symbols), length(entrez),
                  length(unmapped), length(entrez_all) - length(entrez)))
  list(query_name = qname, population = population, gate = gate, direction = direction,
       origin = as.character(spec$origin %||% ""),
       symbols = symbols, entrez = entrez,
       n_symbols = length(symbols), n_mapped = length(entrez_all),
       n_unique_entrez = length(entrez), n_unmapped = length(unmapped),
       unmapped_examples = paste(utils::head(sort(unmapped), 10), collapse = "; "))
}

QUERIES_META <- lapply(QSPEC$sets, build_query)
names(QUERIES_META) <- vapply(QUERIES_META, `[[`, character(1), "query_name")
if (anyDuplicated(names(QUERIES_META)))
  stop("[15_coresh] duplicate query name in coresh.query_signatures.sets.")

QUERIES <- lapply(QUERIES_META, function(q) as.integer(q$entrez))
small   <- vapply(QUERIES, length, integer(1)) < MIN_Q
if (any(small))
  stop("[15_coresh] quer(ies) below min_query_size=", MIN_Q, ": ",
       paste(names(QUERIES)[small], collapse = ", "))

saveRDS(QUERIES, file.path(DIR_OBJECTS, "coresh_hsa_query_entrez.rds"))

# =============================================================================
# 4. Overlap probe -- prove the query actually lands in the compendium
# =============================================================================
# The species/ID failure mode is silent: a mismatched query returns size ~= 0 and
# pctVar noise with no error. Probe ONE chunk first and refuse to sweep if the
# realised per-dataset overlap collapses.

probe <- qs_read(CHUNK_FILES[1])
probe_rows <- rbindlist(lapply(names(QUERIES), function(qn) {
  k <- vapply(probe, function(o) length(stats::na.omit(match(QUERIES[[qn]], o$rownames))), integer(1))
  data.table(query_name = qn, n_query_entrez = length(QUERIES[[qn]]),
             probe_chunk = basename(CHUNK_FILES[1]), n_datasets_probed = length(probe),
             median_overlap = as.numeric(stats::median(k)),
             min_overlap = min(k), max_overlap = max(k),
             median_overlap_frac = as.numeric(stats::median(k)) / length(QUERIES[[qn]]))
}))
print(as.data.frame(probe_rows))
if (any(probe_rows$median_overlap < MIN_Q) || any(probe_rows$median_overlap_frac < 0.05))
  stop("[15_coresh] query/compendium overlap has collapsed (median overlap below ", MIN_Q,
       " genes or 5% of the query). That is the signature of a species or ID-type mismatch, ",
       "NOT a biological null. Refusing to sweep.")
rm(probe); invisible(gc(FALSE))

# =============================================================================
# 5. The sweep (cached) -> tables/coresh_hits.csv
# =============================================================================

load_or_compute <- function(path, fn) {
  if (file.exists(path)) { message("  [cache] ", path); return(readRDS(path)) }
  v <- fn(); saveRDS(v, path); v
}

hits <- load_or_compute(file.path(DIR_OBJECTS, "coresh_hsa_ranked.rds"), function() {
  coresh_batch(queries = QUERIES, chunk_dir = CHUNK_DIR,
               n_cores = N_CORES, pvalues = USE_PVALUES)
})
stopifnot(all(c("query_name", "gse", "gpl", "pctVar", "pval", "size", "rank") %in% names(hits)))
if (!setequal(unique(hits$query_name), names(QUERIES)))
  stop("[15_coresh] cached sweep does not cover the current query matrix. ",
       "Delete ", file.path(DIR_OBJECTS, "coresh_hsa_ranked.rds"), " and re-run.")

hits_dt <- as.data.table(hits)
hits_dt[, `:=`(population = vapply(query_name, function(q) QUERIES_META[[q]]$population, character(1)),
               gate       = vapply(query_name, function(q) QUERIES_META[[q]]$gate,       character(1)))]
setcolorder(hits_dt, c("query_name", "population", "gate", "gse", "gpl", "pctVar", "pval", "size", "rank"))

## NOT a silent cap. The sweep scores EVERY dataset in the compendium, so the full
## ranking is ~265k rows / ~20 MB -- past the point where this repo tracks a table
## (see .gitignore: the tracked surface of 03_results is captions plus compact
## summary/figure-source tables). The exported CSV therefore carries the ranked HEAD
## of each query, the part a reader interprets; the complete ranking stays in the
## checkpoint 03_results/objects/coresh_hsa_ranked.rds and is regenerable from this
## script. Both counts are published in coresh_query_provenance.csv. Raise
## coresh.hits_export_n (or set it to null) to export more.
HITS_EXPORT_N <- CFG$hits_export_n %||% 200L
hits_out <- if (is.null(HITS_EXPORT_N) || is.na(HITS_EXPORT_N)) hits_dt else
  hits_dt[rank <= as.integer(HITS_EXPORT_N)]
data.table::fwrite(hits_out[order(query_name, rank)], file.path(TBL_DIR, "coresh_hits.csv"))
message(sprintf("[15_coresh] sweep: %d rows, %d unique datasets, %d queries; exported top %s per query (%d rows) -> coresh_hits.csv",
                nrow(hits_dt), length(unique(hits_dt$gse)), length(QUERIES),
                as.character(HITS_EXPORT_N %||% "all"), nrow(hits_out)))

# =============================================================================
# 6. Derived co-regulation modules (the CoReSh-to-GSEA bridge)
# =============================================================================

top_hits <- as.data.table(hits)[rank <= TOP_N]
message(sprintf("[15_coresh] deriving modules from the top-%d dataset(s) per query (%d hits).",
                TOP_N, nrow(top_hits)))

gmt_lines <- load_or_compute(file.path(DIR_OBJECTS, "coresh_hsa_gmt_lines.rds"), function() {
  build_coresh_gmt(top_hits = top_hits, queries = QUERIES, chunk_dir = CHUNK_DIR,
                   species = SPECIES, n_top = N_DERIVE,
                   min_size = MIN_SET, max_size = MAX_SET,
                   jaccard_threshold = JACCARD)
})

parse_gmt_lines <- function(lines) {
  if (length(lines) == 0L) return(list())
  out <- lapply(lines, function(ln) {
    tok <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(tok) < 3L) return(NULL)
    unique(tok[-(1:2)][nzchar(tok[-(1:2)])])
  })
  names(out) <- vapply(lines, function(ln) strsplit(ln, "\t", fixed = TRUE)[[1]][1], character(1))
  Filter(function(g) length(g) > 0L, out)
}
DERIVED <- parse_gmt_lines(gmt_lines)
if (length(DERIVED) == 0L)
  stop("[15_coresh] 0 modules survived the size/Jaccard filter -- nothing to score. ",
       "Check query coverage against the compendium.")
saveRDS(DERIVED, file.path(DIR_OBJECTS, "coresh_hsa_derived_sets.rds"))
writeLines(gmt_lines, file.path(TBL_DIR, "coresh_derived_sets.gmt"))

## Set names are "CORESH_<query_name>_<GSE>"; query_name itself contains "_", so
## locate the GSE token rather than splitting positionally.
split_set_name <- function(nm) {
  tok <- strsplit(sub("^CORESH_", "", nm), "_", fixed = TRUE)[[1]]
  gp  <- which(grepl("^GSE", tok))
  gi  <- if (length(gp)) max(gp) else NA_integer_
  list(gse = if (!is.na(gi)) tok[gi] else NA_character_,
       query_name = if (!is.na(gi) && gi > 1L) paste(tok[seq_len(gi - 1L)], collapse = "_") else NA_character_)
}
SET_META <- rbindlist(lapply(names(DERIVED), function(nm) {
  s   <- split_set_name(nm)
  qm  <- QUERIES_META[[s$query_name]]
  hit <- hits_dt[query_name == s$query_name & gse == s$gse]
  seed_syms <- if (is.null(qm)) character(0) else qm$symbols
  data.table(set_name = nm, query_name = s$query_name,
             population = qm$population %||% NA_character_,
             gate = qm$gate %||% NA_character_,
             gse = s$gse, gpl = if (nrow(hit)) hit$gpl[1] else NA_character_,
             hit_rank = if (nrow(hit)) as.integer(hit$rank[1]) else NA_integer_,
             pctVar   = if (nrow(hit)) as.numeric(hit$pctVar[1]) else NA_real_,
             query_genes_in_dataset = if (nrow(hit)) as.integer(hit$size[1]) else NA_integer_,
             n_genes = length(DERIVED[[nm]]),
             n_seed_genes = length(intersect(DERIVED[[nm]], seed_syms)),
             n_new_genes  = length(setdiff(DERIVED[[nm]], seed_syms)))
}))
SET_META[, frac_seed_genes := n_seed_genes / n_genes]
if (any(is.na(SET_META$query_name)))
  stop("[15_coresh] could not parse the seeding query out of one or more module names.")

## Long-form membership table: one row per (module, gene).
members <- rbindlist(lapply(names(DERIVED), function(nm) {
  qm <- QUERIES_META[[split_set_name(nm)$query_name]]
  data.table(set_name = nm, gene = DERIVED[[nm]],
             is_seed_gene = DERIVED[[nm]] %in% qm$symbols)
}))
members <- merge(members, SET_META[, .(set_name, query_name, population, gate, gse, n_genes)],
                 by = "set_name", all.x = TRUE)
setcolorder(members, c("set_name", "query_name", "population", "gate", "gse",
                       "n_genes", "gene", "is_seed_gene"))
data.table::fwrite(members[order(set_name, -is_seed_gene, gene)],
                   file.path(TBL_DIR, "coresh_derived_sets.csv"))
message(sprintf("[15_coresh] modules: %d, sizes %d-%d (median %g) -> coresh_derived_sets.{csv,gmt}",
                nrow(SET_META), min(SET_META$n_genes), max(SET_META$n_genes),
                stats::median(SET_META$n_genes)))

# =============================================================================
# 7. Module-size / ranked-list-overlap distribution
# =============================================================================
# A reader cannot judge a module's enrichment without knowing how many of its
# genes survived into the ranked list being scored -- publish it, per population.

SIZES <- rbindlist(lapply(POPULATIONS, function(p) {
  universe <- names(RANKED[[p]])
  dt <- copy(SET_META)
  dt[, `:=`(scored_against = p,
            n_in_ranked = vapply(set_name, function(nm) length(intersect(DERIVED[[nm]], universe)), integer(1)))]
  dt[, frac_in_ranked := n_in_ranked / n_genes]
  dt[]
}))
data.table::fwrite(SIZES[order(scored_against, set_name)],
                   file.path(TBL_DIR, "coresh_module_sizes.csv"))
## The overlap floor: an empty/NA fgsea result is a silent-failure signature, not
## a null. Refuse to interpret enrichment if the modules do not reach the ranked list.
worst <- SIZES[, .(median_in_ranked = stats::median(n_in_ranked)), by = scored_against]
print(as.data.frame(worst))
if (any(worst$median_in_ranked < GSEA_MIN))
  stop("[15_coresh] median module-to-ranked-list overlap is below gsea_min_size=", GSEA_MIN,
       " for at least one population -- fgsea would return empty results that LOOK like a null. ",
       "Refusing to proceed.")

# =============================================================================
# 8. Pre-ranked GSEA of the modules on the same ranked lists
# =============================================================================
# Engine and parameters are the project's own, read from config -- the same
# clusterProfiler::GSEA(by = "fgsea") call 02_analysis/helpers/fgsea_prerank.R
# makes for the confirmatory scoring stage, so the numbers sit on one scale.

t2g <- data.frame(
  gs_name     = rep(names(DERIVED), lengths(DERIVED)),
  gene_symbol = unlist(DERIVED, use.names = FALSE))

gsea_rows <- list(); gsea_objs <- list()
for (p in POPULATIONS) {
  contrast <- CONTRAST_LABEL[[p]] %||% paste0("SF_vs_PB_", p)
  set.seed(GSEA_SEED)
  g <- clusterProfiler::GSEA(
    geneList = RANKED[[p]], TERM2GENE = t2g, by = "fgsea", exponent = 1, eps = 0,
    minGSSize = GSEA_MIN, maxGSSize = GSEA_MAX, nPermSimple = GSEA_NPERM,
    pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE, verbose = FALSE)
  res <- g@result
  if (nrow(res) == 0L)
    stop("[15_coresh] fgsea returned 0 rows for ", p,
         " despite a verified module-to-ranked-list overlap -- investigate before reading this as a null.")
  gsea_objs[[p]] <- g
  saveRDS(g, file.path(DIR_OBJECTS, sprintf("coresh_hsa_gsea_%s.rds", p)))

  rows <- data.table(
    pathway_id = res$ID, pathway_name = res$ID, database = "CoReSh_derived_hsa",
    nes = res$NES, pvalue = res$pvalue, padj = res$p.adjust,
    set_size = res$setSize, core_enrichment = res$core_enrichment,
    contrast = contrast, population = p,
    direction = ifelse(res$NES >= 0, "Up", "Down"))
  rows <- merge(rows, SET_META[, .(pathway_id = set_name, query_name,
                                   seed_population = population, gate, gse, gpl,
                                   hit_rank, pctVar, n_genes, n_seed_genes, frac_seed_genes)],
                by = "pathway_id", all.x = TRUE)
  rows[, seeded_from_this_population := seed_population == p]
  gsea_rows[[p]] <- rows
  ## No silent caps: a module whose ranked-list overlap fell under gsea_min_size is
  ## dropped by the engine, so name it rather than let the row count quietly shrink.
  dropped <- setdiff(names(DERIVED), rows$pathway_id)
  if (length(dropped))
    message(sprintf("  [%s] %d module(s) dropped below gsea_min_size=%d: %s",
                    p, length(dropped), GSEA_MIN, paste(dropped, collapse = ", ")))
  message(sprintf("  [%s] %d modules scored, %d at padj < %.2f (%d of those seeded from %s itself)",
                  p, nrow(rows), sum(rows$padj < GSEA_FDR, na.rm = TRUE), GSEA_FDR,
                  sum(rows$padj < GSEA_FDR & rows$seeded_from_this_population, na.rm = TRUE), p))
}
GSEA_ALL <- rbindlist(gsea_rows)
setcolorder(GSEA_ALL, c("pathway_id", "database", "population", "contrast", "nes", "pvalue",
                        "padj", "set_size", "direction", "query_name", "seed_population",
                        "seeded_from_this_population", "gate", "gse", "gpl", "hit_rank",
                        "pctVar", "n_genes", "n_seed_genes", "frac_seed_genes"))
data.table::fwrite(GSEA_ALL[order(population, padj)], file.path(TBL_DIR, "coresh_derived_gsea.csv"))
message(sprintf("[15_coresh] derived-set GSEA -> coresh_derived_gsea.csv (%d rows)", nrow(GSEA_ALL)))

## Per-population summary, including the quantity that decides how much of this
## stage is circular: the rank correlation between a module's NES and the fraction
## of that module that is the query that seeded it. A high correlation means the
## enrichment is largely reporting seed content back to itself.
GSEA_SUMMARY <- GSEA_ALL[, .(
  n_modules      = .N,
  n_sig_fdr      = sum(padj < GSEA_FDR, na.rm = TRUE),
  n_sig_up       = sum(padj < GSEA_FDR & nes > 0, na.rm = TRUE),
  n_sig_down     = sum(padj < GSEA_FDR & nes < 0, na.rm = TRUE),
  n_seeded_here  = sum(seeded_from_this_population),
  nes_min        = min(nes), nes_max = max(nes),
  median_set_size      = as.numeric(stats::median(set_size)),
  median_frac_seed     = as.numeric(stats::median(frac_seed_genes)),
  spearman_nes_vs_frac_seed = as.numeric(stats::cor(nes, frac_seed_genes, method = "spearman")),
  top_module     = pathway_id[which.max(nes)], top_nes = max(nes),
  bottom_module  = pathway_id[which.min(nes)], bottom_nes = min(nes)
), by = .(population, contrast)]
data.table::fwrite(GSEA_SUMMARY, file.path(TBL_DIR, "coresh_gsea_summary.csv"))
print(as.data.frame(GSEA_SUMMARY[, .(population, n_sig_up, n_sig_down, nes_min, nes_max,
                                     median_frac_seed, spearman_nes_vs_frac_seed)]))

# =============================================================================
# 9. Running-sum substrate for the top modules (interactive-curve schema)
# =============================================================================
# Schema is FIXED by 03_results/05_scoring/tables/runsum_interactive_treg_WT_heat_up.csv
# and consumed by a downstream notebook:
#   rank, gene, stat, running_es, hit, leading_edge, gene_set, population, contrast
# The running ES is recomputed with the DOSE weighted-KS definition and the fitted
# exponent, so the emitted curve is the one the GSEA object itself describes.
# Mirrors the block in 02_analysis/helpers/fgsea_prerank.R.

n_runsum <- 0L
for (p in POPULATIONS) {
  g   <- gsea_objs[[p]]
  res <- g@result
  contrast <- CONTRAST_LABEL[[p]] %||% paste0("SF_vs_PB_", p)
  exponent <- as.numeric(g@params[["exponent"]] %||% 1)
  gl <- g@geneList; N <- length(gl); gnames <- names(gl)
  top <- res$ID[order(-abs(res$NES))][seq_len(min(RUNSUM_TOP, nrow(res)))]
  for (nm in top) {
    hits_v <- gnames %in% DERIVED[[nm]]
    Nh <- sum(hits_v)
    Phit <- numeric(N); Pmiss <- numeric(N)
    NR <- sum(abs(gl[hits_v])^exponent)
    if (NR > 0) Phit[hits_v] <- (abs(gl[hits_v])^exponent) / NR
    Pmiss[!hits_v] <- 1 / max(N - Nh, 1L)
    running_es <- cumsum(Phit) - cumsum(Pmiss)
    core_str <- res[res$ID == nm, "core_enrichment"]
    core <- if (length(core_str) == 1L && !is.na(core_str) && nzchar(core_str))
      strsplit(core_str, "/")[[1]] else character(0)
    df <- data.frame(
      rank = seq_len(N), gene = gnames, stat = as.numeric(gl),
      running_es = running_es, hit = hits_v,
      leading_edge = hits_v & (gnames %in% core),
      gene_set = nm, population = p, contrast = contrast)
    utils::write.csv(df, file.path(TBL_DIR,
      sprintf("runsum_interactive_coresh_%s_%s.csv", p, nm)), row.names = FALSE)
    n_runsum <- n_runsum + 1L
  }
}
message(sprintf("[15_coresh] running-sum substrate: %d table(s) (top %d modules per population)",
                n_runsum, RUNSUM_TOP))

# =============================================================================
# 10. What the recovered datasets are -- annotation FROM THE CACHE ITSELF
# =============================================================================
# The mouse stage answered this with frozen external web research, because the
# mmu chunks carry only variance structure. The hsa chunks carry more: each
# dataset object holds a `wordMatrix`, the compendium's own centred per-sample
# indicator matrix over the terms that vary most across that dataset's GEO sample
# metadata. That is a checkable, in-cache descriptor of what a dataset contrasts,
# so it is used here instead of external research -- and named for what it is.
# Correlating each term against the query direction says which metadata term
# tracks the axis the query defines in that dataset. Descriptive only: it never
# enters a statistic and cannot change any pctVar or NES.

## GSE -> chunk index. .index_chunks() (from the skill's bridge) memoises in
## options() within a session, but building it means reading all 89 chunks, so it
## is checkpointed here too -- otherwise a re-run with a warm sweep cache still
## pays the full scan.
gse_index <- load_or_compute(file.path(DIR_OBJECTS, "coresh_hsa_gse_index.rds"),
                             function() .index_chunks(CHUNK_DIR))
want <- unique(SET_META[, .(gse, query_name)])
ann_rows <- list()
for (ck in unique(gse_index[gse %in% want$gse]$chunk)) {
  ch <- qs_read(ck)
  ids <- vapply(ch, function(o) o$gseId, character(1))
  for (i in which(ids %in% want$gse)) {
    o <- ch[[i]]
    wm <- o$wordMatrix
    for (qn in unique(want[gse == o$gseId]$query_name)) {
      qi <- stats::na.omit(match(QUERIES[[qn]], o$rownames))
      prof <- colSums((o$E1024 / 1024)[qi, , drop = FALSE])
      prof <- prof / sqrt(sum(prof^2))
      r <- if (is.null(wm) || is.null(dim(wm)) || nrow(wm) == 0L) numeric(0) else
        suppressWarnings(apply(wm, 1L, function(w)
          if (stats::sd(w) == 0) NA_real_ else stats::cor(w, prof)))
      ord <- if (length(r)) order(abs(r), decreasing = TRUE, na.last = TRUE) else integer(0)
      trm <- names(r)[ord]; rr <- unname(r[ord])
      keep <- seq_len(min(N_TERMS, length(trm)))
      ann_rows[[length(ann_rows) + 1L]] <- data.table(
        gse = o$gseId, ann_gpl = o$gplId, recovering_query = qn,
        organism = "Homo sapiens", n_samples = as.integer(o$nsamples),
        n_genes_profiled = length(o$rownames),
        metadata_terms = paste(trm[keep], collapse = " / "),
        terms_aligned_with_query = paste(utils::head(trm[keep][rr[keep] > 0 & !is.na(rr[keep])], 5), collapse = " / "),
        terms_opposed_to_query   = paste(utils::head(trm[keep][rr[keep] < 0 & !is.na(rr[keep])], 5), collapse = " / "),
        max_abs_term_r = if (length(rr)) max(abs(rr), na.rm = TRUE) else NA_real_)
    }
  }
  rm(ch); invisible(gc(FALSE))
}
ANN_RAW <- rbindlist(ann_rows)
## The compendium is keyed by (GSE, GPL): a series measured on two platforms is two
## objects, so a bare GSE join would double the annotation. Prefer the platform the
## ranked hit actually scored on, keep the multiplicity visible, and fall back to the
## first record only if that platform is somehow absent.
ANN_RAW[, n_platform_records_in_compendium := .N, by = .(gse, recovering_query)]
key <- SET_META[, .(set_name, gse, gpl, recovering_query = query_name, population, gate,
                    hit_rank, pctVar, n_genes, n_seed_genes, frac_seed_genes)]
ANN <- merge(key, ANN_RAW, by.x = c("gse", "gpl", "recovering_query"),
             by.y = c("gse", "ann_gpl", "recovering_query"), all.x = TRUE)
if (any(is.na(ANN$n_samples))) {
  fb <- ANN_RAW[, .SD[1L], by = .(gse, recovering_query)]
  need <- ANN[is.na(n_samples), .(set_name, gse, recovering_query)]
  message(sprintf("  [annotation] %d set(s) had no record on the ranked platform; using the first (%s).",
                  nrow(need), paste(need$set_name, collapse = ", ")))
  ANN <- rbind(ANN[!is.na(n_samples)],
               merge(key[set_name %in% need$set_name],
                     fb[, !c("ann_gpl"), with = FALSE], by = c("gse", "recovering_query")),
               fill = TRUE)
}
if (nrow(ANN) != nrow(SET_META))
  stop("[15_coresh] annotation join produced ", nrow(ANN), " rows for ", nrow(SET_META), " modules.")
ANN[, `:=`(annotation_source = "CoReSh preprocessed chunk wordMatrix (per-sample GEO metadata term matrix)",
           annotation_basis  = "Pearson r between each term's centred per-sample indicator and the normalised query profile",
           cache_snapshot    = CACHE_TAG, cache_synapse_id = CACHE_SYN, cache_downloaded_at = CACHE_DATE)]
data.table::fwrite(ANN[order(population, gate, hit_rank)],
                   file.path(TBL_DIR, "coresh_derived_annotation.csv"))
message(sprintf("[15_coresh] dataset annotation -> coresh_derived_annotation.csv (%d rows)", nrow(ANN)))

# =============================================================================
# 11. Query provenance -- what was asked, and what the ID seam cost
# =============================================================================

PROV <- rbindlist(lapply(QUERIES_META, function(q) {
  pr <- probe_rows[query_name == q$query_name]
  data.table(
    query_name = q$query_name, population = q$population, gate = q$gate,
    direction = q$direction, origin = q$origin,
    contrast = CONTRAST_LABEL[[q$population]] %||% NA_character_,
    source_ranked_list = fill_tmpl(RNK_TMPL, q$population),
    n_symbols_in_ranked_list = length(RANKED[[q$population]]),
    gate_rule = if (identical(q$gate, "fdr_logfc"))
      sprintf("padj < %g AND log2FoldChange >= %g", DE_FDR, DE_LOGFC)
      else sprintf("padj < %g AND log2FoldChange > 0", DE_FDR),
    gate_source_table = fill_tmpl(DE_TMPL, q$population),
    n_query_symbols = q$n_symbols,
    n_mapped_to_entrez = q$n_mapped,
    n_unique_entrez = q$n_unique_entrez,
    n_dropped_unmapped = q$n_unmapped,
    frac_dropped_unmapped = q$n_unmapped / q$n_symbols,
    n_collapsed_duplicate_entrez = q$n_mapped - q$n_unique_entrez,
    unmapped_examples = q$unmapped_examples,
    id_map_package = paste0("org.Hs.eg.db ", as.character(utils::packageVersion("org.Hs.eg.db"))),
    species = SPECIES, chunk_dir = CHUNK_DIR, n_chunks = length(CHUNK_FILES),
    n_datasets_searched = length(unique(hits_dt[query_name == q$query_name]$gse)),
    n_hits_exported_to_csv = nrow(hits_out[query_name == q$query_name]),
    full_ranking_checkpoint = file.path(DIR_OBJECTS, "coresh_hsa_ranked.rds"),
    cache_synapse_id = CACHE_SYN, cache_snapshot = CACHE_TAG, cache_downloaded_at = CACHE_DATE,
    probe_median_overlap = if (nrow(pr)) pr$median_overlap else NA_real_,
    probe_median_overlap_frac = if (nrow(pr)) pr$median_overlap_frac else NA_real_,
    search_pvalues = USE_PVALUES, top_n_hits = TOP_N, n_derive = N_DERIVE,
    derived_min_size = MIN_SET, derived_max_size = MAX_SET, jaccard_threshold = JACCARD,
    gsea_engine = "clusterProfiler::GSEA(by='fgsea')",
    gsea_min_size = GSEA_MIN, gsea_max_size = GSEA_MAX,
    gsea_nperm = GSEA_NPERM, gsea_seed = GSEA_SEED,
    fgsea_version = as.character(utils::packageVersion("fgsea")),
    script = SCRIPT_PATH)
}))
data.table::fwrite(PROV, file.path(TBL_DIR, "coresh_query_provenance.csv"))
message("[15_coresh] provenance -> coresh_query_provenance.csv")

# =============================================================================
# 12. Structural asserts
# =============================================================================

expected <- c("coresh_query_provenance.csv", "coresh_hits.csv", "coresh_derived_sets.csv",
              "coresh_derived_sets.gmt", "coresh_derived_gsea.csv", "coresh_gsea_summary.csv",
              "coresh_derived_annotation.csv", "coresh_module_sizes.csv")
missing <- expected[!file.exists(file.path(TBL_DIR, expected))]
if (length(missing)) stop("[15_coresh] missing output(s): ", paste(missing, collapse = ", "))
stopifnot(n_runsum > 0L,
          nrow(GSEA_ALL) <= length(DERIVED) * length(POPULATIONS),
          nrow(GSEA_ALL) > 0L,
          all(!is.na(GSEA_ALL$padj)))

message(sprintf(
  "[15_coresh] COMPLETE -- %d queries over %d public human datasets, %d modules, %d GSEA rows. Exploratory tier.",
  length(QUERIES), length(unique(hits_dt$gse)), length(DERIVED), nrow(GSEA_ALL)))
